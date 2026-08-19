#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "../../third_party/get_max_cpus/c_codes/get_max_cpus.hpp"
#include "read_and_write_data.tpp"
#include "sz3_hohs_decompress.tpp"
#include "validate_params.hpp"

namespace mrcz
{
    namespace detail
    {
        template<typename T,typename Backend>
        void decompress_field(param_class<T> &params,const string &field,Backend backend)
        {
            double br=0,cr=0,max_err=0,mse=0,psnr=0;
            ifstream compressed_file(params.data_path+"_"+field+"."+params.method,ios::binary);
            if (!compressed_file)
                throw runtime_error("无法打开压缩文件: "+params.data_path+"_"+field+"."+params.method);
            compressed_file.seekg(1,ios::beg);
            unsigned char is_constant=0;
            compressed_file.read(reinterpret_cast<char*>(&is_constant),1);
            if (is_constant)
            {
                compressed_file.read(reinterpret_cast<char*>(&params.data_range_min[field]),sizeof(T));
                params.data_range_max[field]=params.data_range_min[field];
                params.decompressed_data.assign(params.element_num_to_process*params.header.element_size,params.data_range_min[field]);
                br=static_cast<double>(compressed_file.tellg())*8/(params.element_num_to_process*params.header.element_size);
                cr=(params.header.word_size*8)/br;
            }
            else
            {
                uint64_t partition_num=0;
                compressed_file.read(reinterpret_cast<char*>(&partition_num),sizeof(partition_num));
                params.partitions.resize(partition_num);
                for (uint64_t i=0;i<partition_num;i++)
                {
                    auto &partition=params.partitions[i];
                    partition.partition_id=i;
                    partition.begin_index=(i*params.element_num_to_process)/partition_num;
                    partition.end_index=((i+1)*params.element_num_to_process)/partition_num;
                    partition.element_num=partition.end_index-partition.begin_index;
                }
                compressed_file.read(reinterpret_cast<char*>(&params.data_range_min[field]),sizeof(T));
                compressed_file.read(reinterpret_cast<char*>(&params.data_range_max[field]),sizeof(T));
                for (auto &partition:params.partitions)
                {
                    uint64_t compressed_size=0;
                    compressed_file.read(reinterpret_cast<char*>(&compressed_size),sizeof(compressed_size));
                    partition.zstd_bs.resize(compressed_size);
                }
                for (auto &partition:params.partitions)
                    compressed_file.read(reinterpret_cast<char*>(partition.zstd_bs.data()),partition.zstd_bs.size());
                br=static_cast<double>(compressed_file.tellg())*8/(params.element_num_to_process*params.header.element_size);
                cr=(params.header.word_size*8)/br;
                compressed_file.close();

                T data_range=params.data_range_max[field]-params.data_range_min[field];
                if (params.eb_mode=="ABS")
                {
                    params.abs_eb[field]=params.input_eb;
                    params.rel_eb[field]=params.input_eb/data_range;
                }
                else
                {
                    params.rel_eb[field]=params.input_eb;
                    params.abs_eb[field]=params.input_eb*data_range;
                }
                printf("params.method= %s\n",params.method.c_str());
                backend(params,field);
                for (auto &partition:params.partitions)
                {
                    partition.zstd_bs.clear();
                    partition.zstd_bs.shrink_to_fit();
                }
            }

            if (params.fields_to_process.size()!=params.all_fields_in_file.size() or
                params.header.file_element_num!=params.element_num_to_process)
                printf("没有选择全部的field，或者没有选择全部的element进行解压，不生成解压文件\n");
            else if constexpr (is_same_v<T,double>)
            {
                if (params.header.word_size==4)
                {
                    param_class<float> float_params(params);
                    float_params.decompressed_data.assign(params.decompressed_data.begin(),params.decompressed_data.end());
                    write_data(float_params,field);
                }
                else
                    write_data(params,field);
            }
            else
                write_data(params,field);

            ifstream original_file(params.data_path,ios::binary);
            if (original_file)
            {
                original_file.close();
                params.data.resize(params.element_num_to_process*params.header.element_size);
                if constexpr (is_same_v<T,double>)
                {
                    if (params.header.word_size==4)
                    {
                        param_class<float> float_params(params);
                        float_params.data.resize(params.data.size());
                        ::read_data(float_params,field);
                        params.data.assign(float_params.data.begin(),float_params.data.end());
                    }
                    else
                        ::read_data(params,field);
                }
                else
                    ::read_data(params,field);

                #pragma omp parallel for reduction(max:max_err) reduction(+:mse) schedule(static) num_threads(params.max_cpus)
                for (uint64_t i=0;i<params.data.size();i++)
                {
                    T error=fabs(params.data[i]-params.decompressed_data[i]);
                    max_err=max(max_err,static_cast<double>(error));
                    mse+=error*error;
                }
                if (!is_constant and max_err>params.abs_eb[field]*1.1)
                {
                    uint64_t max_err_index=0;
                    for (uint64_t i=0;i<params.data.size();i++)
                        if (fabs(params.data[i]-params.decompressed_data[i])==max_err)
                        {
                            max_err_index=i;
                            break;
                        }
                    printf("Warning: index=%lu, original=%.17e, decompressed=%.17e, max_err=%e, abs_eb=%e\n",
                        max_err_index,static_cast<double>(params.data[max_err_index]),
                        static_cast<double>(params.decompressed_data[max_err_index]),max_err,params.abs_eb[field]);
                }
                mse/=params.data.size();
                double data_range=params.data_range_max[field]-params.data_range_min[field];
                psnr=(mse==0) ? INFINITY : 10*log10(data_range*data_range/mse);
            }
            else
                printf("原始数据不存在，无法计算压缩质量\n");

            if (is_constant)
            {
                params.abs_eb[field]=0;
                params.rel_eb[field]=0;
            }
            printf("\033[91mField= %s, abs_eb= %.15f, rel_eb= %.15f, BR= %.15f, CR= %.15f, max_err= %.15f (%f%%), MSE= %.15f, PSNR= %.15f dB\033[0m\n",
                field.c_str(),params.abs_eb[field],params.rel_eb[field],br,cr,max_err,
                params.abs_eb[field]==0 ? 0 : max_err/params.abs_eb[field]*100,mse,psnr);
            params.data.clear();
            params.decompressed_data.clear();
        }

        template<typename T,typename Backend>
        void decompress_all_fields_with_backend(param_class<T> &params,Backend backend)
        {
            for (auto &field:params.fields_to_process)
            {
                printf("Field= %s\n",field.c_str());
                ifstream compressed_file(params.data_path+"_"+field+"."+params.method,ios::binary);
                if (!compressed_file)
                    throw runtime_error("无法打开压缩文件: "+params.data_path+"_"+field+"."+params.method);
                uint8_t compressed_word_size=0;
                compressed_file.read(reinterpret_cast<char*>(&compressed_word_size),1);
                if constexpr (is_same_v<T,float>)
                {
                    if (compressed_word_size==8)
                    {
                        param_class<double> double_params(params);
                        decompress_field(double_params,field,backend);
                    }
                    else
                        decompress_field(params,field,backend);
                }
                else
                    decompress_field(params,field,backend);
            }
        }

        template<typename T,typename Backend>
        void decompress_file_with_backend(param_class<T> &params,Backend backend)
        {
            validate_common_params(params);
            if (params.method.empty())
                params.method="mrcz-3d-gll";
            params.all_fields_to_process=params.fields_to_process.empty();
            params.read_header("decompress");
            params.default_init();
            params.save_header_and_element_map("decompress");
            params.max_cpus=get_max_cpus(params.max_cpus);
            if (params.header.word_size==8 and is_same_v<T,float>)
            {
                param_class<double> double_params(params);
                decompress_all_fields_with_backend(double_params,backend);
            }
            else
                decompress_all_fields_with_backend(params,backend);
        }
    }

    namespace core
    {
        template<typename T>
        void decompress(param_class<T> &params)
        {
            if (params.method.empty())
                params.method="mrcz-3d-gll";
            if (params.method!="mrcz-3d-gll" and params.method!="mrcz-3d-gll-experimental" and params.method!="sz3-3d-gll")
                throw runtime_error("mrcz::core::decompress不支持方法: "+params.method);
            auto backend=[]<typename U>(param_class<U> &backend_params,string field)
            {
                if (backend_params.method=="sz3-3d-gll")
                    sz3_2_decompress(backend_params,field);
                else
                    backend_params.decompress_mrcz(field);
            };
            detail::decompress_file_with_backend(params,backend);
        }
    }
}
