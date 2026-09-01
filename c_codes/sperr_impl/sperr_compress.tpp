#pragma once

#include "../baseline_bitstream.hpp"
#include "../third_party/sperr/include/SPERR_C_API.h"

template<typename T>
vector<unsigned char> sperr_compress_field(mrcz::param_class<T> &params,mrcz::partition_class &partition,string method,T *data,T abs_eb)
{
    auto &header=params.header;
    size_t n0=header.c_element_shape[0];
    size_t n1=header.c_element_shape[1];
    size_t n2=header.c_element_shape[2];
    int is_float=is_same<T,float>::value?1:0;
    void *compressed=nullptr;
    size_t compressed_size=0;
    int status;
    if (method=="SPERR-2D")
        status=C_API::sperr_comp_2d(data,is_float,header.element_size,partition.element_num,3,abs_eb,0,&compressed,&compressed_size);
    else if (method=="SPERR-3D-fused" or method=="SPERR-3D-per-element")
    {
        size_t element_num=method=="SPERR-3D-per-element"?1:partition.element_num;
        status=C_API::sperr_comp_3d(data,is_float,n2,n1,element_num*n0,n2,n1,n0,3,abs_eb,1,&compressed,&compressed_size);
    }
    else
        throw runtime_error("未知的SPERR布局: "+method);
    if (status!=0)
        throw runtime_error("SPERR压缩失败，状态码="+to_string(status));
    vector<unsigned char> output(reinterpret_cast<unsigned char*>(compressed),reinterpret_cast<unsigned char*>(compressed)+compressed_size);
    free(compressed);
    return output;
}

template<typename T>
void sperr_compress(mrcz::param_class<T> &params,string field)
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
        if (params.method=="SPERR-3D-per-element")
        {
            partition.zstd_bs.clear();
            for (uint64_t element_id=0;element_id<partition.element_num;element_id++)
            {
                vector<unsigned char> compressed=sperr_compress_field(params,partition,params.method,data.data()+begin+element_id*header.element_size,abs_eb);
                mrcz::baseline::append_uint64(partition.zstd_bs,compressed.size());
                partition.zstd_bs.insert(partition.zstd_bs.end(),compressed.begin(),compressed.end());
            }
        }
        else
            partition.zstd_bs=sperr_compress_field(params,partition,params.method,data.data()+begin,abs_eb);
    }
    double end_time=omp_get_wtime();
    printf("Compression time= %f seconds\n",end_time-start_time);
    printf("\033[91mField= %s, compression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),params.element_num_to_process*header.element_size*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
}
