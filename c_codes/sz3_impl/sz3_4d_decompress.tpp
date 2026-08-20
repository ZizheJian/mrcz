#pragma once

#include "../mrcz_core/mrcz.hpp"
#include "SZ3/api/sz.hpp"

template<typename T>
void sz3_4d_decompress(mrcz::param_class<T> &params,string field)
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
        if (partition.zstd_bs.empty())
            throw runtime_error("SZ3-4D-uniform partition bitstream为空");

        uint64_t begin=partition.begin_index*header.element_size;
        SZ3::Config config;
        T *output=decompressed_data.data()+begin;
        SZ_decompress<T>(config,reinterpret_cast<const char*>(partition.zstd_bs.data()),
            partition.zstd_bs.size(),output);
        if (config.N!=4 or config.dims[0]!=partition.element_num or
            config.dims[1]!=header.element_shape[2] or config.dims[2]!=header.element_shape[1] or
            config.dims[3]!=header.element_shape[0])
            throw runtime_error("SZ3-4D-uniform partition尺寸与bitstream不一致");
    }
    double end_time=omp_get_wtime();
    printf("Decompression time= %f seconds\n",end_time-start_time);
    printf("\033[91mField= %s, decompression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),
        params.element_num_to_process*header.element_size*header.word_size/
        ((end_time-start_time)*1024*1024*params.max_cpus));
}
