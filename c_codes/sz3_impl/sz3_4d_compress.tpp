#pragma once

#include "../mrcz_core/mrcz.hpp"
#include "SZ3/api/sz.hpp"

template<typename T>
void sz3_4d_compress(mrcz::param_class<T> &params,string field)
{
    auto &header=params.header;
    auto &partitions=params.partitions;
    auto &data=params.data;
    T abs_eb=params.abs_eb[field];
    uint32_t n0=header.c_element_shape[0];
    uint32_t n1=header.c_element_shape[1];
    uint32_t n2=header.c_element_shape[2];

    double start_time=omp_get_wtime();
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
    for (uint64_t partition_id=0;partition_id<partitions.size();partition_id++)
    {
        auto &partition=partitions[partition_id];
        uint64_t begin=partition.begin_index*header.element_size;

        // SZ3的最后一个维度变化最快，因此[element,n2,n1,n0]与MRCZ的内存顺序一致。
        SZ3::Config config(partition.element_num,n2,n1,n0);
        config.errorBoundMode=SZ3::EB_ABS;
        config.absErrorBound=abs_eb;
        size_t compressed_size=0;
        char *compressed_data=SZ_compress<T>(config,data.data()+begin,compressed_size);
        partition.zstd_bs.assign(compressed_data,compressed_data+compressed_size);
        delete[] compressed_data;
    }
    double end_time=omp_get_wtime();
    printf("Compression time= %f seconds\n",end_time-start_time);
    printf("\033[91mField= %s, compression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),
        params.element_num_to_process*header.element_size*header.word_size/
        ((end_time-start_time)*1024*1024*params.max_cpus));
}
