#pragma once

#include <set>
#include "../mrcz_core/mrcz_param.hpp"
#include "../mrcz_core/huffman.hpp"
#include "../mrcz_core/quantize.tpp"
#include "../third_party/zstd/lib/zstd.h"


template<typename T>
void sz2_1_compress(mrcz::param_class<T> &params,string field)
{
    auto &header=params.header;
    auto &blocks=params.blocks;
    auto &data=params.data;
    auto &x0=params.x[0];
    auto &x1=params.x[1];
    auto &x2=params.x[2];

    uint32_t n0=header.element_shape[0];
    uint32_t n1=header.element_shape[1];
    uint32_t n2=header.element_shape[2];
    uint32_t n01=n0*n1;
    uint32_t n02=n0*n2;
    uint32_t n12=n1*n2;
    uint32_t n012=n0*n1*n2;
    T abs_eb=params.abs_eb[field];
    T neg_2_eb=1/(2*abs_eb);
    T coeff_eb=0.1*abs_eb;
    T neg_2_coeff_eb=1/(2*coeff_eb);
    int64_t *err_qb=new int64_t[params.element_num_to_process*(n012+4)];
    T sum_x00=0,sum_x11=0,sum_x22=0;
    for (uint32_t i=0;i<x0.size();i++)
        sum_x00+=x0[i]*x0[i];
    for (uint32_t i=0;i<x1.size();i++)
        sum_x11+=x1[i]*x1[i];
    for (uint32_t i=0;i<x2.size();i++)
        sum_x22+=x2[i]*x2[i];
    double start_time=omp_get_wtime();
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
    for (uint64_t block_id=0;block_id<blocks.size();block_id++)
    {
        auto &block=blocks[block_id];
        uint64_t begin_index=block.begin_index;
        uint64_t end_index=block.end_index;
        int64_t *a_address=err_qb+begin_index*(n012+4);
        int64_t *b_address=a_address+end_index-begin_index;
        int64_t *c_address=b_address+end_index-begin_index;
        int64_t *d_address=c_address+end_index-begin_index;
        int64_t *block_err_qb=d_address+end_index-begin_index;

        vector<T> sum_element_0(n12,0),sum_element_1(n02,0),sum_element_01(n2,0),sum_element_02(n1,0),sum_element_12(n0,0);
        T a,b,c,d;
        for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
        {
            T *element=data.data()+element_id*n012;
            sum_element_0.assign(n12,0);
            sum_element_1.assign(n02,0);
            sum_element_01.assign(n2,0);
            sum_element_02.assign(n1,0);
            sum_element_12.assign(n0,0);
            a=b=c=d=0;
            for (uint32_t i0=0;i0<n0;i0++)
                for (uint32_t i1=0;i1<n1;i1++)
                    for (uint32_t i2=0;i2<n2;i2++)
                    {
                        sum_element_0[i1*n2+i2]+=element[i0*n12+i1*n2+i2];
                        sum_element_1[i0*n2+i2]+=element[i0*n12+i1*n2+i2];
                    }
            for (uint32_t i1=0;i1<n1;i1++)
                for (uint32_t i2=0;i2<n2;i2++)
                {
                    sum_element_01[i2]+=sum_element_0[i1*n2+i2];
                    sum_element_02[i1]+=sum_element_0[i1*n2+i2];
                }
            for (uint32_t i0=0;i0<n0;i0++)
                for (uint32_t i2=0;i2<n2;i2++)
                    sum_element_12[i0]+=sum_element_1[i0*n2+i2];
            if (sum_x00!=0)
            {
                for (uint32_t i0=0;i0<n0;i0++)
                    a+=x0[i0]*sum_element_12[i0];
                a/=(n12*sum_x00);
            }
            else
                a=0;
            if (sum_x11!=0)
            {
                for (uint32_t i1=0;i1<n1;i1++)
                    b+=x1[i1]*sum_element_02[i1];
                b/=(n02*sum_x11);
            }
            else
                b=0;
            if (sum_x22!=0)
            {
                for (uint32_t i2=0;i2<n2;i2++)
                    c+=x2[i2]*sum_element_01[i2];
                c/=(n01*sum_x22);
            }
            else
                c=0;
            for (uint32_t i0=0;i0<n0;i0++)
                d+=sum_element_12[i0];
            d/=n012;
            params.quantize_and_overwrite(a,a_address[element_id-begin_index],coeff_eb,neg_2_coeff_eb);
            params.quantize_and_overwrite(b,b_address[element_id-begin_index],coeff_eb,neg_2_coeff_eb);
            params.quantize_and_overwrite(c,c_address[element_id-begin_index],coeff_eb,neg_2_coeff_eb);
            params.quantize_and_overwrite(d,d_address[element_id-begin_index],coeff_eb,neg_2_coeff_eb);

            vector<T> new_element(n012,0);
            for (uint32_t i0=0;i0<n0;i0++)
                for (uint32_t i1=0;i1<n1;i1++)
                    for (uint32_t i2=0;i2<n2;i2++)
                        new_element[i0*n12+i1*n2+i2]=a*x0[i0]+b*x1[i1]+c*x2[i2]+d;
            params.quantize_err_array(new_element.data(),element,block_err_qb+(element_id-begin_index)*n012,n012,abs_eb,neg_2_eb);
        }
        // // 把两种qb无损压缩到一起
        vector<unsigned char> huffman_bs_vec((n012+4)*(end_index-begin_index)*sizeof(int64_t));
        unsigned char *huffman_bs=huffman_bs_vec.data();
        mrcz::huffman_encoder_class<int64_t> err_encoder;
        err_encoder.preprocess_encode(a_address,(n012+4)*(end_index-begin_index));
        err_encoder.save(huffman_bs);
        err_encoder.encode(a_address,(n012+4)*(end_index-begin_index),huffman_bs);
        uint64_t huffman_bs_len=huffman_bs-huffman_bs_vec.data();
        uint64_t max_zstd_bs_len=ZSTD_compressBound(huffman_bs_len);
        block.zstd_bs.resize(max_zstd_bs_len);
        uint64_t zstd_bs_len=ZSTD_compress(block.zstd_bs.data(),max_zstd_bs_len,huffman_bs_vec.data(),huffman_bs_len,3);
        block.zstd_bs.resize(zstd_bs_len);
    }
    double end_time=omp_get_wtime();
    delete[] err_qb;
    printf("Compression time= %f seconds\n",end_time-start_time);
    printf("\033[91mBandwidth per cpu= %f MB/s\n\033[0m",params.element_num_to_process*n012*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
}