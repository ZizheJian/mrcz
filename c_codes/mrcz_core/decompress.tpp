#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <get_max_cpus.hpp>
#include "read_and_write_data.tpp"
#include "fdwz_decompress.tpp"
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
            unsigned char constant_flag=0;
            compressed_file.read(reinterpret_cast<char*>(&constant_flag),1);
            if (constant_flag==1)
            {
                T constant_value=0;
                compressed_file.read(reinterpret_cast<char*>(&constant_value),sizeof(T));
                params.data_range_min[field]=params.data_range_max[field]=constant_value;
                params.decompressed_data.assign(params.element_num_to_process*params.header.element_size,constant_value);
                br=static_cast<double>(compressed_file.tellg())*8/(params.element_num_to_process*params.header.element_size);
                cr=(params.header.word_size*8)/br;
            }
            else if (constant_flag==0)
            {
                uint64_t partition_num=0;
                compressed_file.read(reinterpret_cast<char*>(&partition_num),sizeof(partition_num));
                if (partition_num==0 or partition_num>params.element_num_to_process)
                    throw runtime_error("bitstream中的partition数量无效: "+to_string(partition_num));
                if (params.partition_num!=0 and params.partition_num!=partition_num)
                    throw runtime_error("命令行指定的partition数量与bitstream不一致: requested="+to_string(params.partition_num)+", stored="+to_string(partition_num));
                params.partition_num=partition_num;
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
                printf("Max CPUs= %d, partition num= %lu\n",params.max_cpus,params.partition_num);

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
            else
                throw runtime_error("不支持的constant flag: "+to_string(constant_flag));

            if (params.fields_to_process.size()!=params.all_fields_in_file.size() or
                params.header.file_element_num!=params.element_num_to_process)
                printf("没有选择全部的field，或者没有选择全部的element进行解压，不生成解压文件\n");
            else if constexpr (is_same_v<T,double>)
            {
                if (params.header.word_size==4)
                {
                    param_class<float> float_params(params,"decompressed_data");
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
                uint64_t value_num=params.element_num_to_process*params.header.element_size;
                if constexpr (is_same_v<T,double>)
                {
                    if (params.header.word_size==4)
                    {
                        param_class<float> float_params(params);
                        float_params.data.resize(value_num);
                        ::read_data(float_params,field);
                        params.data.assign(float_params.data.begin(),float_params.data.end());
                    }
                    else
                    {
                        params.data.resize(value_num);
                        ::read_data(params,field);
                    }
                }
                else
                {
                    params.data.resize(value_num);
                    ::read_data(params,field);
                }

                if (constant_flag==1)
                {
                    auto range=minmax_element(params.data.begin(),params.data.end());
                    params.data_range_min[field]=*range.first;
                    params.data_range_max[field]=*range.second;
                    T data_range=*range.second-*range.first;
                    if (params.eb_mode=="ABS")
                    {
                        params.abs_eb[field]=params.input_eb;
                        params.rel_eb[field]=data_range==0 ? 0 : params.input_eb/data_range;
                    }
                    else
                    {
                        params.rel_eb[field]=params.input_eb;
                        params.abs_eb[field]=params.input_eb*data_range;
                    }
                }

                uint64_t zero_error_num=0,error_ge_50_num=0,error_ge_90_num=0,error_ge_99_num=0,error_ge_999_num=0,error_ge_9999_num=0,error_ge_100_num=0,error_gt_10001_num=0;
                constexpr uint32_t error_histogram_bin_num=10000;
                vector<vector<uint64_t>> error_histograms;
                if(params.if_analyze_error_distribution)
                    error_histograms.assign(params.max_cpus,vector<uint64_t>(error_histogram_bin_num+1,0));
                #pragma omp parallel for reduction(max:max_err) reduction(+:mse,zero_error_num,error_ge_50_num,error_ge_90_num,error_ge_99_num,error_ge_999_num,error_ge_9999_num,error_ge_100_num,error_gt_10001_num) schedule(static) num_threads(params.max_cpus)
                for (uint64_t i=0;i<params.data.size();i++)
                {
                    T error=fabs(params.data[i]-params.decompressed_data[i]);
                    max_err=max(max_err,static_cast<double>(error));
                    mse+=error*error;
                    if(params.if_analyze_error_distribution)
                    {
                        double normalized_error=static_cast<double>(error)/params.abs_eb[field];
                        zero_error_num+=error==0;
                        error_ge_50_num+=normalized_error>=0.5;
                        error_ge_90_num+=normalized_error>=0.9;
                        error_ge_99_num+=normalized_error>=0.99;
                        error_ge_999_num+=normalized_error>=0.999;
                        error_ge_9999_num+=normalized_error>=0.9999;
                        error_ge_100_num+=normalized_error>=1;
                        error_gt_10001_num+=normalized_error>1.0001;
                        uint32_t bin=normalized_error>=1?error_histogram_bin_num:static_cast<uint32_t>(normalized_error*error_histogram_bin_num);
                        error_histograms[omp_get_thread_num()][bin]++;
                    }
                }
                if(params.if_analyze_error_distribution)
                {
                    vector<uint64_t> error_histogram(error_histogram_bin_num+1,0);
                    for(auto &thread_histogram:error_histograms)
                        for(uint32_t bin=0;bin<=error_histogram_bin_num;bin++)
                            error_histogram[bin]+=thread_histogram[bin];
                    auto quantile_upper_bound=[&](double quantile)
                    {
                        uint64_t target=max(uint64_t(1),static_cast<uint64_t>(ceil(quantile*params.data.size()))),cumulative=0;
                        for(uint32_t bin=0;bin<=error_histogram_bin_num;bin++)
                        {
                            cumulative+=error_histogram[bin];
                            if(cumulative>=target)
                                return bin==error_histogram_bin_num?max_err/params.abs_eb[field]:static_cast<double>(bin+1)/error_histogram_bin_num;
                        }
                        return max_err/params.abs_eb[field];
                    };
                    double value_num=params.data.size();
                    printf("Error distribution: field=%s, value_num=%lu, zero=%lu (%.9e), p50_upper=%.7f, p90_upper=%.7f, p99_upper=%.7f, p999_upper=%.7f, p9999_upper=%.7f, p99999_upper=%.7f, ge_0.5=%lu (%.9e), ge_0.9=%lu (%.9e), ge_0.99=%lu (%.9e), ge_0.999=%lu (%.9e), ge_0.9999=%lu (%.9e), ge_1=%lu (%.9e), gt_1.0001=%lu (%.9e)\n",
                        field.c_str(),params.data.size(),zero_error_num,zero_error_num/value_num,quantile_upper_bound(0.5),quantile_upper_bound(0.9),quantile_upper_bound(0.99),quantile_upper_bound(0.999),quantile_upper_bound(0.9999),quantile_upper_bound(0.99999),
                        error_ge_50_num,error_ge_50_num/value_num,error_ge_90_num,error_ge_90_num/value_num,error_ge_99_num,error_ge_99_num/value_num,error_ge_999_num,error_ge_999_num/value_num,error_ge_9999_num,error_ge_9999_num/value_num,error_ge_100_num,error_ge_100_num/value_num,error_gt_10001_num,error_gt_10001_num/value_num);
                }
                if (max_err>params.abs_eb[field]*1.1)
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
                params.method="MRCZ";
            params.all_fields_to_process=params.fields_to_process.empty();
            params.read_header("decompress");
            params.default_init();
            params.save_header_and_element_map("decompress");
            params.max_cpus=get_max_cpus(static_cast<int>(min<uint64_t>(params.max_cpus,params.element_num_to_process)));
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
                params.method="MRCZ";
            if (params.method=="MRCZ")
            {
                params.grid_assumption="GLL";
                params.transform="DChT";
                params.coeff_eb_scale=0.05;
                params.if_coefficient_transpose=true;
                params.if_tiled_coefficient_transpose=true;
            }
            else if (params.method=="FDWZ" or params.method=="FDWZ-ablation")
                params.grid_assumption="GLL";
            auto backend=[]<typename U>(param_class<U> &backend_params,string field)
            {
                if (backend_params.method=="MRCZ")
                {
                    if (backend_params.if_measure_phase_time)
                        backend_params.template decompress_mrcz<true,false>(field);
                    else
                        backend_params.template decompress_mrcz<false,false>(field);
                }
                else if (backend_params.method=="MRCZ-ablation")
                    backend_params.template decompress_mrcz<false,true>(field);
                else if (backend_params.method=="FDWZ")
                {
                    if(backend_params.if_measure_phase_time)
                        backend_params.template decompress_fdwz<true,false>(field);
                    else
                        backend_params.template decompress_fdwz<false,false>(field);
                }
                else if (backend_params.method=="FDWZ-ablation")
                    backend_params.template decompress_fdwz<false,true>(field);
                else
                    throw runtime_error("未知的core解压方法: "+backend_params.method);
            };
            detail::decompress_file_with_backend(params,backend);
        }
    }
}
