#pragma once

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "../../third_party/get_max_cpus/c_codes/get_max_cpus.hpp"
#include "read_and_write_data.tpp"
#include "sz3_hohs_compress.tpp"
#include "validate_params.hpp"

namespace mrcz
{
    namespace detail
    {
        template<typename T,typename Backend>
        void compress_field_and_write(param_class<T> &params,const string &field,Backend backend)
        {
            bool is_constant=params.data_range_min[field]==params.data_range_max[field];
            if (is_constant)
                printf("Warning: field %s has constant value %f\n",field.c_str(),params.data_range_min[field]);
            else
                backend(params,field);

            double compressed_bytes=0;
            string compressed_path=params.data_path+"_"+field+"."+params.method;
            ofstream compressed_file(compressed_path,ios::binary);
            if (!compressed_file)
                throw runtime_error("无法创建压缩文件: "+compressed_path);

            uint8_t word_size_during_compression=sizeof(T);
            compressed_file.write(reinterpret_cast<const char*>(&word_size_during_compression),1);
            if (is_constant)
            {
                unsigned char constant_flag=1;
                compressed_file.write(reinterpret_cast<const char*>(&constant_flag),1);
                compressed_file.write(reinterpret_cast<const char*>(&params.data_range_min[field]),sizeof(T));
                compressed_bytes=2+sizeof(T);
            }
            else
            {
                unsigned char constant_flag=0;
                compressed_file.write(reinterpret_cast<const char*>(&constant_flag),1);
                uint64_t partition_num=params.partitions.size();
                compressed_file.write(reinterpret_cast<const char*>(&partition_num),sizeof(partition_num));
                compressed_file.write(reinterpret_cast<const char*>(&params.data_range_min[field]),sizeof(T));
                compressed_file.write(reinterpret_cast<const char*>(&params.data_range_max[field]),sizeof(T));
                compressed_bytes=2+sizeof(partition_num)+2*sizeof(T);
                for (auto &partition:params.partitions)
                {
                    uint64_t compressed_size=partition.zstd_bs.size();
                    compressed_file.write(reinterpret_cast<const char*>(&compressed_size),sizeof(compressed_size));
                    compressed_bytes+=sizeof(compressed_size);
                }
                for (auto &partition:params.partitions)
                {
                    compressed_file.write(reinterpret_cast<const char*>(partition.zstd_bs.data()),partition.zstd_bs.size());
                    compressed_bytes+=partition.zstd_bs.size();
                    partition.zstd_bs.clear();
                    partition.zstd_bs.shrink_to_fit();
                }
            }
            compressed_file.close();

            double bits_per_value=compressed_bytes*8/(params.element_num_to_process*params.header.element_size);
            double compression_ratio=(8*params.header.word_size)/bits_per_value;
            printf("\033[91mBR= %f bits/value, CR= %f\n\033[0m",bits_per_value,compression_ratio);
        }

        template<typename T,typename Backend>
        void compress_all_fields_with_backend(param_class<T> &params,Backend backend)
        {
            for (auto &field:params.fields_to_process)
            {
                printf("Field= %s\n",field.c_str());
                params.data.resize(params.element_num_to_process*params.header.element_size);
                ::read_data(params,field);

                T global_min=params.data[0];
                T global_max=params.data[0];
                #pragma omp parallel for schedule(static) num_threads(params.max_cpus) reduction(min:global_min) reduction(max:global_max)
                for (uint64_t partition_id=0;partition_id<params.partitions.size();partition_id++)
                {
                    auto &partition=params.partitions[partition_id];
                    uint64_t begin_address=partition.begin_index*params.header.element_size;
                    uint64_t end_address=partition.end_index*params.header.element_size;
                    for (uint64_t i=begin_address;i<end_address;i++)
                    {
                        global_min=min(global_min,params.data[i]);
                        global_max=max(global_max,params.data[i]);
                    }
                }
                params.data_range_min[field]=global_min;
                params.data_range_max[field]=global_max;
                T data_range=global_max-global_min;
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
                printf("Data range: [%f, %f], abs_eb=%e, rel_eb=%e\n",global_min,global_max,params.abs_eb[field],params.rel_eb[field]);

                if constexpr (is_same_v<T,float>)
                {
                    if (params.rel_eb[field]<1e-4)
                    {
                        param_class<double> double_params(params);
                        if (double_params.eb_mode=="ABS")
                        {
                            double_params.abs_eb[field]=double_params.input_eb;
                            double_params.rel_eb[field]=double_params.input_eb/
                                (double_params.data_range_max[field]-double_params.data_range_min[field]);
                        }
                        else
                        {
                            double_params.rel_eb[field]=double_params.input_eb;
                            double_params.abs_eb[field]=double_params.input_eb*
                                (double_params.data_range_max[field]-double_params.data_range_min[field]);
                        }
                        double_params.data.assign(params.data.begin(),params.data.end());
                        params.data.clear();
                        params.data.shrink_to_fit();
                        compress_field_and_write(double_params,field,backend);
                    }
                    else
                        compress_field_and_write(params,field,backend);
                }
                else
                    compress_field_and_write(params,field,backend);

                params.data.clear();
                params.data.shrink_to_fit();
            }
        }

        template<typename T,typename Backend>
        void compress_file_with_backend(param_class<T> &params,Backend backend)
        {
            validate_common_params(params);
            params.all_fields_to_process=params.fields_to_process.empty();
            params.read_header("compress");
            params.default_init();
            params.save_header_and_element_map("compress");
            params.max_cpus=get_max_cpus(params.max_cpus);

            params.partitions.resize(params.max_cpus);
            for (int i=0;i<params.max_cpus;i++)
            {
                auto &partition=params.partitions[i];
                partition.partition_id=i;
                partition.begin_index=(i*params.element_num_to_process)/params.max_cpus;
                partition.end_index=((i+1)*params.element_num_to_process)/params.max_cpus;
                partition.element_num=partition.end_index-partition.begin_index;
            }

            if (params.header.word_size==8 and is_same_v<T,float>)
            {
                param_class<double> double_params(params);
                compress_all_fields_with_backend(double_params,backend);
            }
            else
                compress_all_fields_with_backend(params,backend);
        }
    }

    namespace core
    {
        template<typename T>
        void compress(param_class<T> &params)
        {
            if (params.method.empty())
                params.method="mrcz";
            if (params.method!="mrcz" and params.method!="mrcz-experimental" and params.method!="sz3-hohs")
                throw runtime_error("mrcz::core::compress不支持方法: "+params.method);

            auto backend=[]<typename U>(param_class<U> &backend_params,string field)
            {
                if (backend_params.method=="sz3-hohs")
                    sz3_2_compress(backend_params,field);
                else
                    backend_params.compress_mrcz(field);
            };
            detail::compress_file_with_backend(params,backend);
        }
    }
}
