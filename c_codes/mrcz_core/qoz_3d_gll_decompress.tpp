#pragma once

#include "huffman.hpp"
#include "mrcz.hpp"
#include "qoz_3d_gll_prediction.tpp"
#include "../third_party/zstd/lib/zstd.h"

template<typename T>
void qoz_3d_gll_decompress(mrcz::param_class<T> &params,string field)
{
    auto &header=params.header;
    auto &blocks=params.partitions;
    auto &decompressed_data=params.decompressed_data;
    uint32_t n0=header.c_element_shape[0],n1=header.c_element_shape[1],n2=header.c_element_shape[2];
    uint32_t n12=n1*n2,n012=n0*n12,corner_num=1;
    for (uint32_t dim=0;dim<3;dim++)
        if (header.c_element_shape[dim]>1)
            corner_num*=2;
    T abs_eb=params.abs_eb[field];
    int64_t *corner_qb=new int64_t[corner_num*params.element_num_to_process];
    int64_t *err_qb=new int64_t[params.element_num_to_process*n012];
    int64_t *no_corner_err_qb=new int64_t[params.element_num_to_process*(n012-corner_num)];
    decompressed_data.resize(params.element_num_to_process*n012);
    mrcz::detail::qoz_3d_gll_prediction_plan<T> prediction_plan(params.x);

    double start_time=omp_get_wtime();
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
    for (uint64_t block_id=0;block_id<blocks.size();block_id++)
    {
        auto &block=blocks[block_id];
        uint64_t begin_index=block.begin_index,end_index=block.end_index;
        uint64_t huffman_bs_len=ZSTD_getFrameContentSize(block.zstd_bs.data(),block.zstd_bs.size());
        vector<unsigned char> huffman_bs_vec(huffman_bs_len);
        if (ZSTD_decompress(huffman_bs_vec.data(),huffman_bs_len,block.zstd_bs.data(),block.zstd_bs.size())!=huffman_bs_len)
            throw runtime_error("ZSTD_decompress error in qoz_3d_gll_decompress()");
        const unsigned char *huffman_bs=huffman_bs_vec.data();
        mrcz::huffman_encoder_class<int64_t> corner_encoder;
        corner_encoder.load(huffman_bs);
        corner_encoder.decode(huffman_bs,corner_qb+corner_num*begin_index,corner_num*(end_index-begin_index));
        mrcz::huffman_encoder_class<int64_t> err_encoder;
        err_encoder.load(huffman_bs);
        err_encoder.decode(huffman_bs,no_corner_err_qb+(n012-corner_num)*begin_index,(n012-corner_num)*(end_index-begin_index));
        params.reconstruct_corner_number(corner_qb,begin_index,end_index,n0,n1,n2,n12,n012,abs_eb);
        for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
        {
            int64_t *element_err_qb=err_qb+element_id*n012;
            int64_t *element_no_corner_err_qb=no_corner_err_qb+element_id*(n012-corner_num);
            params.convert_no_corner_err_eb_to_err_eb(element_no_corner_err_qb,element_err_qb,n0,n1,n2,n12,n012);
            T *element=decompressed_data.data()+element_id*n012;
            auto apply_prediction=[&](T prediction,T &truth,int64_t &qb)
            {
                truth=prediction-2*abs_eb*qb;
            };
            mrcz::detail::qoz_3d_gll_process_element(element,element_err_qb,prediction_plan,n1,n2,apply_prediction);
        }
    }
    double end_time=omp_get_wtime();
    delete[] corner_qb;
    delete[] err_qb;
    delete[] no_corner_err_qb;
    printf("Decompression time= %f seconds\n",end_time-start_time);
    printf("\033[91mField= %s, decompression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),params.element_num_to_process*n012*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
}
