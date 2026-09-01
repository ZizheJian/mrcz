#pragma once

#include "../baseline_bitstream.hpp"
#include "SZ3/api/sz.hpp"

template<typename T>
SZ3::Config make_sz3_api_config(mrcz::param_class<T> &params,mrcz::partition_class &partition,string method)
{
    auto &header=params.header;
    uint64_t n0=header.c_element_shape[0];
    uint64_t n1=header.c_element_shape[1];
    uint64_t n2=header.c_element_shape[2];
    if (method=="SZ3-1D")
        return SZ3::Config(partition.element_num*header.element_size);
    if (method=="SZ3-2D")
        return SZ3::Config(partition.element_num,header.element_size);
    if (method=="SZ3-3D-fused")
        return SZ3::Config(partition.element_num*n0,n1,n2);
    if (method=="SZ3-4D")
        return SZ3::Config(partition.element_num,n0,n1,n2);
    if (method=="SZ3-3D-per-element")
        return SZ3::Config(n0,n1,n2);
    throw runtime_error("未知的SZ3布局: "+method);
}

template<typename T>
void sz3_compress(mrcz::param_class<T> &params,string field)
{
    auto &header=params.header;
    auto &partitions=params.partitions;
    auto &data=params.data;
    T abs_eb=params.abs_eb[field];
    double start_time=omp_get_wtime();
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
    for (uint64_t partition_id=0;partition_id<partitions.size();partition_id++)
    {
        auto &partition=partitions[partition_id];
        uint64_t begin=partition.begin_index*header.element_size;
        if (params.method=="SZ3-3D-per-element")
        {
            partition.zstd_bs.clear();
            for (uint64_t element_id=0;element_id<partition.element_num;element_id++)
            {
                SZ3::Config config=make_sz3_api_config(params,partition,params.method);
                config.errorBoundMode=SZ3::EB_ABS;
                config.absErrorBound=abs_eb;
                size_t compressed_size=0;
                char *compressed_data=SZ_compress<T>(config,data.data()+begin+element_id*header.element_size,compressed_size);
                mrcz::baseline::append_uint64(partition.zstd_bs,compressed_size);
                partition.zstd_bs.insert(partition.zstd_bs.end(),reinterpret_cast<unsigned char*>(compressed_data),reinterpret_cast<unsigned char*>(compressed_data)+compressed_size);
                delete[] compressed_data;
            }
        }
        else
        {
            SZ3::Config config=make_sz3_api_config(params,partition,params.method);
            config.errorBoundMode=SZ3::EB_ABS;
            config.absErrorBound=abs_eb;
            size_t compressed_size=0;
            char *compressed_data=SZ_compress<T>(config,data.data()+begin,compressed_size);
            partition.zstd_bs.assign(compressed_data,compressed_data+compressed_size);
            delete[] compressed_data;
        }
    }
    double end_time=omp_get_wtime();
    printf("Compression time= %f seconds\n",end_time-start_time);
    printf("\033[91mField= %s, compression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),params.element_num_to_process*header.element_size*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
}
