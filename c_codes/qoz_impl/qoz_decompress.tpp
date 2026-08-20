#pragma once

#include "../mrcz_core/mrcz.hpp"
#include "qoz_bridge.hpp"

template<typename T>
void qoz_decompress(mrcz::param_class<T> &params,string field)
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
            throw runtime_error("QoZ-4D-uniform partition bitstream为空");
        uint64_t begin=partition.begin_index*header.element_size;
        qoz_decompress_4d(partition.zstd_bs,decompressed_data.data()+begin);
    }
    double end_time=omp_get_wtime();
    printf("Decompression time= %f seconds\n",end_time-start_time);
    printf("\033[91mField= %s, decompression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),
        params.element_num_to_process*header.element_size*header.word_size/
        ((end_time-start_time)*1024*1024*params.max_cpus));
}
