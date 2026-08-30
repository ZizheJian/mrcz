#pragma once

#include "../baseline_api_layout.hpp"
#include "sz.h"

struct sz2_api_dims
{
    size_t r5=0,r4=0,r3=0,r2=0,r1=0;
};

template<typename T>
sz2_api_dims make_sz2_api_dims(mrcz::param_class<T> &params,mrcz::partition_class &partition,mrcz::baseline::api_layout layout)
{
    auto &header=params.header;
    size_t n0=header.c_element_shape[0];
    size_t n1=header.c_element_shape[1];
    size_t n2=header.c_element_shape[2];
    if (layout==mrcz::baseline::api_layout::flattened_1d)
        return {.r1=partition.element_num*header.element_size};
    if (layout==mrcz::baseline::api_layout::elements_2d)
        return {.r2=partition.element_num,.r1=header.element_size};
    if (layout==mrcz::baseline::api_layout::fused_3d)
        return {.r3=partition.element_num*n0,.r2=n1,.r1=n2};
    if (layout==mrcz::baseline::api_layout::elements_4d)
        return {.r4=partition.element_num,.r3=n0,.r2=n1,.r1=n2};
    return {.r3=n0,.r2=n1,.r1=n2};
}

template<typename T>
vector<unsigned char> sz2_compress_field(mrcz::param_class<T> &params,mrcz::partition_class &partition,mrcz::baseline::api_layout layout,T *data,T abs_eb)
{
    sz2_api_dims dims=make_sz2_api_dims(params,partition,layout);
    size_t compressed_size=0;
    int data_type=is_same<T,float>::value?SZ_FLOAT:SZ_DOUBLE;
    unsigned char *compressed=SZ_compress_args(data_type,data,&compressed_size,ABS,abs_eb,0,0,dims.r5,dims.r4,dims.r3,dims.r2,dims.r1);
    if (compressed==nullptr)
        throw runtime_error("SZ2压缩失败");
    vector<unsigned char> output(compressed,compressed+compressed_size);
    free(compressed);
    return output;
}

template<typename T>
void sz2_api_compress(mrcz::param_class<T> &params,string field)
{
    auto &header=params.header;
    auto &partitions=params.partitions;
    auto &data=params.data;
    T abs_eb=params.abs_eb[field];
    auto layout=mrcz::baseline::parse_api_layout(params.method,"SZ2");
    if (SZ_Init(nullptr)!=SZ_SCES)
        throw runtime_error("SZ2初始化失败");
    double start_time=omp_get_wtime();
    // XXX SZ2官方API在一次进程内共享全局状态，并发压缩多个partition会产生损坏的bitstream
    for (uint64_t partition_id=0;partition_id<partitions.size();partition_id++)
    {
        auto &partition=partitions[partition_id];
        uint64_t begin=partition.begin_index*header.element_size;
        if (layout==mrcz::baseline::api_layout::per_element_3d)
        {
            partition.zstd_bs.clear();
            mrcz::baseline::append_uint64(partition.zstd_bs,partition.element_num);
            for (uint64_t element_id=0;element_id<partition.element_num;element_id++)
            {
                vector<unsigned char> compressed=sz2_compress_field(params,partition,layout,data.data()+begin+element_id*header.element_size,abs_eb);
                mrcz::baseline::append_uint64(partition.zstd_bs,compressed.size());
                partition.zstd_bs.insert(partition.zstd_bs.end(),compressed.begin(),compressed.end());
            }
        }
        else
            partition.zstd_bs=sz2_compress_field(params,partition,layout,data.data()+begin,abs_eb);
    }
    double end_time=omp_get_wtime();
    SZ_Finalize();
    printf("Compression time= %f seconds\n",end_time-start_time);
    printf("\033[91mField= %s, compression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),params.element_num_to_process*header.element_size*header.word_size/((end_time-start_time)*1024*1024));
}
