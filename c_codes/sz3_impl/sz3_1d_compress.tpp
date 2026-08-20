#pragma once

#include "../mrcz_core/mrcz.hpp"
#include "SZ3/api/sz.hpp"

template<typename T>
void sz3_1d_compress(mrcz::param_class<T> &params,string field)
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
        uint64_t value_num=partition.element_num*header.element_size;

        SZ3::Config config(value_num);
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
