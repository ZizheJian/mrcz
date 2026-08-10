#pragma once

#include "../mrcz_core/mrcz_param.hpp"
#include "../third_party/zfp/include/zfp.h"

template<typename T>
void zfp_decompress(mrcz::param_class<T> &params,string field)
{
	auto &header=params.header;
    auto &blocks=params.partitions;
    auto &decompressed_data=params.decompressed_data;

    uint32_t n0=header.element_shape[0];
    uint32_t n1=header.element_shape[1];
    uint32_t n2=header.element_shape[2];
    T abs_eb=params.abs_eb[field];
    zfp_type type=std::is_same<T,float>::value ? zfp_type_float : zfp_type_double;
    decompressed_data.resize(params.element_num_to_process*header.element_size);
    double start_time=omp_get_wtime();
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
    for (uint64_t block_id=0;block_id<blocks.size();block_id++)
    {
        auto &block=blocks[block_id];
        uint64_t begin_index=block.begin_index;
        uint64_t end_index=block.end_index;
        
        zfp_stream *zfp=zfp_stream_open(NULL);
        zfp_stream_set_accuracy(zfp,abs_eb);
        bitstream *input_bs=stream_open(block.zstd_bs.data(),block.zstd_bs.size());
        zfp_stream_set_bit_stream(zfp,input_bs);
        zfp_field *field=zfp_field_4d(decompressed_data.data()+begin_index*header.element_size,type,n2,n1,n0,end_index-begin_index);
        zfp_decompress(zfp,field);
        stream_close(input_bs);
        zfp_field_free(field);
        zfp_stream_close(zfp);
    }
    double end_time=omp_get_wtime();
    printf("Decompression time= %f seconds\n",end_time-start_time);
    printf("\033[91mField= %s, decompression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),params.element_num_to_process*header.element_size*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
}