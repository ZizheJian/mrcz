#pragma once


#include "../mrcz_core/mrcz.hpp"
#include "../third_party/sperr/include/SPERR_C_API.h"

template<typename T>
void sperr_decompress(mrcz::param_class<T> &params,string field)
{
	auto &header=params.header;
    auto &blocks=params.partitions;
    auto &decompressed_data=params.decompressed_data;

    int is_float=is_same<T,float>::value?1:0;
    decompressed_data.resize(params.element_num_to_process*header.element_size);
    double start_time=omp_get_wtime();
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
    for (uint64_t block_id=0;block_id<blocks.size();block_id++)
    {
        auto &block=blocks[block_id];
        uint64_t begin_index=block.begin_index;
        uint64_t end_index=block.end_index;

        void *decompressed_buffer=nullptr;
        size_t dimx,dimy,dimz;
        C_API::sperr_decomp_3d(block.zstd_bs.data(),block.zstd_bs.size(),is_float,0,(size_t*)&dimx,(size_t*)&dimy,(size_t*)&dimz,&decompressed_buffer);
        memcpy(decompressed_data.data()+begin_index*header.element_size,decompressed_buffer,(end_index-begin_index)*header.element_size*sizeof(T));
        free(decompressed_buffer);
    }
    double end_time=omp_get_wtime();
    printf("Decompression time= %f seconds\n",end_time-start_time);
    printf("\033[91mField= %s, decompression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),params.element_num_to_process*header.element_size*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
}
