#pragma once

#include "../mrcz_core/mrcz.hpp"
#include "SZ3/api/sz.hpp"

template<typename T>
void sz3_1d_decompress(mrcz::param_class<T> &params,string field)
{
    auto &header=params.header;
    auto &partitions=params.partitions;
    auto &decompressed_data=params.decompressed_data;
    decompressed_data.resize(params.element_num_to_process*header.element_size);

    double start_time=omp_get_wtime();
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
    for (uint64_t partition_id=0;partition_id<partitions.size();partition_id++)
    {
        auto &partition=partitions[partition_id];
        uint64_t begin=partition.begin_index*header.element_size;
        if (partition.zstd_bs.empty())
            throw runtime_error("SZ3-1D-uniform partition bitstream为空");

        SZ3::Config config;
        T *output=decompressed_data.data()+begin;
        SZ_decompress<T>(config,reinterpret_cast<const char*>(partition.zstd_bs.data()),
            partition.zstd_bs.size(),output);
        if (config.N!=1 or config.num!=partition.element_num*header.element_size)
            throw runtime_error("SZ3-1D-uniform partition尺寸与bitstream不一致");
    }
    double end_time=omp_get_wtime();
    printf("Decompression time= %f seconds\n",end_time-start_time);
    printf("\033[91mField= %s, decompression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),
        params.element_num_to_process*header.element_size*header.word_size/
        ((end_time-start_time)*1024*1024*params.max_cpus));
}
