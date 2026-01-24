#pragma once

#include <set>
#include "../jzzz_core/jzzz_param.hpp"
#include "../jzzz_core/huffman.hpp"
#include "../jzzz_core/quantize.tpp"
#include "../third_party/zstd/lib/zstd.h"


template<typename T>
void sz2_1_compress(jzzz::param_class<T> &params,string field)
{
    auto &header=params.header;
    auto &blocks=params.blocks;
    auto &data=params.data;
    auto &x=params.x;

    uint32_t n=header.element_shape[0];
    uint32_t n2=n*n;
    uint32_t n3=n2*n;
    T abs_eb=params.abs_eb[field];
    T neg_2_eb=1/(2*abs_eb);
    T coeff_eb=0.1*abs_eb;
    T neg_2_coeff_eb=1/(2*coeff_eb);
    int64_t *err_qb=new int64_t[params.element_num_to_process*(header.element_size+4)];
    T sum_x2=0;
    for (uint32_t i=0;i<x.size();i++)
        sum_x2+=x[i]*x[i];
    double start_time=omp_get_wtime();
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
    for (uint64_t block_id=0;block_id<blocks.size();block_id++)
    {
        auto &block=blocks[block_id];
        uint64_t begin_index=block.begin_index;
        uint64_t end_index=block.end_index;
        int64_t *a_address=err_qb+begin_index*(header.element_size+4);
        int64_t *b_address=a_address+end_index-begin_index;
        int64_t *c_address=b_address+end_index-begin_index;
        int64_t *d_address=c_address+end_index-begin_index;
        int64_t *block_err_qb=d_address+end_index-begin_index;

        vector<T> sum_element_x(100,0),sum_element_y(100,0),sum_element_xy(10,0),sum_element_yz(10,0),sum_element_xz(10,0);
        T a,b,c,d;
        for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
        {
            T *element=data.data()+element_id*header.element_size;
            sum_element_x.assign(100,0);
            sum_element_y.assign(100,0);
            sum_element_xy.assign(10,0);
            sum_element_yz.assign(10,0);
            sum_element_xz.assign(10,0);
            a=b=c=d=0;
            for (uint32_t i0=0;i0<n;i0++)
                for (uint32_t i1=0;i1<n;i1++)
                    for (uint32_t i2=0;i2<n;i2++)
                    {
                        sum_element_x[i1*n+i2]+=element[i0*n2+i1*n+i2];
                        sum_element_y[i0*n+i2]+=element[i0*n2+i1*n+i2];
                    }
            for (uint32_t i0=0;i0<n;i0++)
                for (uint32_t i1=0;i1<n;i1++)
                {
                    sum_element_xy[i1]+=sum_element_x[i0*n+i1];
                    sum_element_yz[i0]+=sum_element_y[i0*n+i1];
                    sum_element_xz[i0]+=sum_element_x[i0*n+i1];
                }
            for (uint32_t i0=0;i0<n;i0++)
                a+=x[i0]*sum_element_yz[i0];
            a/=(n*n*sum_x2);
            for (uint32_t i0=0;i0<n;i0++)
                b+=x[i0]*sum_element_xz[i0];
            b/=(n*n*sum_x2);
            for (uint32_t i0=0;i0<n;i0++)
                c+=x[i0]*sum_element_xy[i0];
            c/=(n*n*sum_x2);
            for (uint32_t i0=0;i0<n;i0++)
                d+=sum_element_xy[i0];
            d/=(n3);
            params.quantize_and_overwrite(a,a_address[element_id-begin_index],coeff_eb,neg_2_coeff_eb);
            params.quantize_and_overwrite(b,b_address[element_id-begin_index],coeff_eb,neg_2_coeff_eb);
            params.quantize_and_overwrite(c,c_address[element_id-begin_index],coeff_eb,neg_2_coeff_eb);
            params.quantize_and_overwrite(d,d_address[element_id-begin_index],coeff_eb,neg_2_coeff_eb);

            vector<T> new_element(n3,0);
            for (uint32_t i0=0;i0<n;i0++)
                for (uint32_t i1=0;i1<n;i1++)
                    for (uint32_t i2=0;i2<n;i2++)
                        new_element[i0*n2+i1*n+i2]=a*x[i0]+b*x[i1]+c*x[i2]+d;
            params.quantize_err_array(new_element.data(),element,block_err_qb+(element_id-begin_index)*header.element_size,header.element_size,abs_eb,neg_2_eb);
        }
        // // 把两种qb无损压缩到一起
        vector<unsigned char> huffman_bs_vec((header.element_size+4)*(end_index-begin_index)*sizeof(int64_t));
        unsigned char *huffman_bs=huffman_bs_vec.data();
        jzzz::huffman_encoder_class<int64_t> err_encoder;
        err_encoder.preprocess_encode(a_address,(header.element_size+4)*(end_index-begin_index));
        err_encoder.save(huffman_bs);
        err_encoder.encode(a_address,(header.element_size+4)*(end_index-begin_index),huffman_bs);
        uint64_t huffman_bs_len=huffman_bs-huffman_bs_vec.data();
        uint64_t max_zstd_bs_len=ZSTD_compressBound(huffman_bs_len);
        block.zstd_bs.resize(max_zstd_bs_len);
        uint64_t zstd_bs_len=ZSTD_compress(block.zstd_bs.data(),max_zstd_bs_len,huffman_bs_vec.data(),huffman_bs_len,3);
        block.zstd_bs.resize(zstd_bs_len);
    }
    double end_time=omp_get_wtime();
    delete[] err_qb;
    printf("Compression time= %f seconds\n",end_time-start_time);
    printf("\033[91mBandwidth per cpu= %f MB/s\n\033[0m",params.element_num_to_process*header.element_size*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
}