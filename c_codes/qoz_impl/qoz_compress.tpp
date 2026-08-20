#pragma once

#include "../mrcz_core/mrcz.hpp"
#include "qoz_bridge.hpp"

template<typename T>
void qoz_compress(mrcz::param_class<T> &params,string field)
{
    auto &header=params.header;
    auto &partitions=params.partitions;
    auto &data=params.data;
    T abs_eb=params.abs_eb[field];
    uint32_t n0=header.element_shape[0];
    uint32_t n1=header.element_shape[1];
    uint32_t n2=header.element_shape[2];

    double start_time=omp_get_wtime();
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
    for (uint64_t partition_id=0;partition_id<partitions.size();partition_id++)
    {
        auto &partition=partitions[partition_id];
        uint64_t begin=partition.begin_index*header.element_size;
        qoz_compress_4d(data.data()+begin,n0,n1,n2,partition.element_num,abs_eb,partition.zstd_bs);
    }
    double end_time=omp_get_wtime();
    printf("Compression time= %f seconds\n",end_time-start_time);
    printf("\033[91mField= %s, compression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),
        params.element_num_to_process*header.element_size*header.word_size/
        ((end_time-start_time)*1024*1024*params.max_cpus));
}
