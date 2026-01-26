#pragma once

#include "../jzzz_core/jzzz_param.hpp"
#include "../third_party/sperr/include/SPERR_C_API.h"

template<typename T>
void sperr_compress(jzzz::param_class<T> &params,string field)
{
    auto &header=params.header;
    auto &blocks=params.blocks;
    auto &data=params.data;

    uint32_t n=header.element_shape[0];
    T abs_eb=params.abs_eb[field];
    int is_float=std::is_same<T,float>::value?1:0;
    double start_time=omp_get_wtime();
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
    for (uint64_t block_id=0;block_id<blocks.size();block_id++)
    {
        auto &block=blocks[block_id];
        uint64_t begin_index=block.begin_index;
        uint64_t end_index=block.end_index;

        void *temp_zstd_bs=nullptr;
        size_t block_zstd_bs_len=0;
        C_API::sperr_comp_3d((const void*)(data.data()+begin_index*header.element_size),
            is_float,n,n,(end_index-begin_index)*n,n,n,n,3,abs_eb,0,&temp_zstd_bs,&block_zstd_bs_len);
        block.zstd_bs.resize(block_zstd_bs_len);
        memcpy(block.zstd_bs.data(),temp_zstd_bs,block_zstd_bs_len);
        free(temp_zstd_bs);
    }
    double end_time=omp_get_wtime();
    printf("Compression time= %f seconds\n",end_time-start_time);
    printf("\033[91mBandwidth per cpu= %f MB/s\n\033[0m",params.element_num_to_process*header.element_size*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
}