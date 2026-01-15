#ifndef __COMPRESS_CPP__
#define __COMPRESS_CPP__

#include <omp.h>
#include <immintrin.h>
#include "jzzz_param.hpp"
#include "huffman.cpp"
#include <cmath>
#include "../third_party/zstd/lib/zstd.h"

namespace jzzz
{
    void param_class::compress()
    {
        const char *slurm_cpus=getenv("SLURM_CPUS_PER_TASK");
        int available_cpus=0;
        if (slurm_cpus!=nullptr)
            available_cpus=atoi(slurm_cpus);
        else
            throw runtime_error("没有获取到每个task可用的cpu数量");
        printf("Available cpus per task= %d\n",available_cpus);
        uint32_t n=header.element_shape[0];
        param_eb=eb*param_eb_scale;
        float zero_freq_bound=param_eb*zero_f_qb_bound;
        float neg_2_param_eb=1/(2*param_eb);
        float neg_2_eb=1/(2*eb);
        int64_t *corner_qb=new int64_t[8*element_num_to_process];
        int64_t *f_qb=new int64_t[element_num_to_process*header.element_size];
        int64_t *err_qb=new int64_t[element_num_to_process*header.element_size];
        zstd_bs.resize(blocks.size());
        double start_time=omp_get_wtime();
        #pragma omp parallel for schedule(static) num_threads(available_cpus)
        for (uint64_t fake_block_id=0;fake_block_id<blocks_to_process.size();fake_block_id++) //这个block_id是假的，指的是在需要处理的blocks中的编号，而不是header文件中的编号
        {
            uint64_t block_id=blocks_to_process[fake_block_id];
            printf("CPU %d compressing block %ld\n",omp_get_thread_num(),block_id);
            uint64_t begin_index=blocks[block_id].begin_index;
            uint64_t end_index=blocks[block_id].end_index;
            // 先量化位于每个element位于角上的点，覆盖原数据
            // 之后可能要改成压缩率更高的方法
            vector<float> corner_numbers(8*(end_index-begin_index));
            for (uint64_t i=begin_index;i<end_index;i++)
            {
                float *element=data.data()+i*header.element_size;
                int64_t *element_corner_qb=corner_qb+8*i;
                quantize_and_overwrite_number(element[0],element_corner_qb[0],eb,neg_2_eb);
                quantize_and_overwrite_number(element[n-1],element_corner_qb[1],eb,neg_2_eb);
                quantize_and_overwrite_number(element[(n-1)*n],element_corner_qb[2],eb,neg_2_eb);
                quantize_and_overwrite_number(element[n*n-1],element_corner_qb[3],eb,neg_2_eb);
                quantize_and_overwrite_number(element[(n-1)*n*n],element_corner_qb[4],eb,neg_2_eb);
                quantize_and_overwrite_number(element[(n-1)*n*n+n-1],element_corner_qb[5],eb,neg_2_eb);
                quantize_and_overwrite_number(element[(n-1)*n*n+(n-1)*n],element_corner_qb[6],eb,neg_2_eb);
                quantize_and_overwrite_number(element[n*n*n-1],element_corner_qb[7],eb,neg_2_eb);
            }
            // 然后量化剩余的数据
            for (uint64_t i=begin_index;i<end_index;i++)
            {
                float *element=data.data()+i*header.element_size;
                int64_t *element_corner_qb=corner_qb+8*i;
                int64_t *element_f_qb=f_qb+i*header.element_size;
                int64_t *element_err_qb=err_qb+i*header.element_size;
                __m256 vec_corner_numbers=_mm256_set_ps(element[n*n*n-1],element[(n-1)*n*n+(n-1)*n],element[(n-1)*n*n+n-1],element[(n-1)*n*n],
                                                        element[n*n-1],element[(n-1)*n],element[n-1],element[0]);
                for (uint32_t i0=0;i0<n;i0++)
                    for (uint32_t i1=0;i1<n;i1++)
                        for (uint32_t i2=0;i2<n;i2++)
                        {
                            uint32_t pos=i0*n*n+i1*n+i2;
                            uint32_t bias=8*pos;
                            __m256 vec_element_base_weight=_mm256_loadu_ps(&element_base_weight[bias]);
                            __m256 vec_mult=_mm256_mul_ps(vec_corner_numbers,vec_element_base_weight);
                            __m256 sum1=_mm256_hadd_ps(vec_mult,vec_mult);
                            __m256 sum2=_mm256_hadd_ps(sum1,sum1);
                            __m128 low=_mm256_castps256_ps128(sum2);
                            __m128 high=_mm256_extractf128_ps(sum2,1);
                            __m128 sum3=_mm_add_ps(low,high);
                            float total_weight=_mm_cvtss_f32(sum3);
                            element[pos]-=total_weight;
                        }
                vector<float> temp_element1(n*n*n,0),temp_element2(n*n*n,0);
                for (uint32_t i=0;i<n;i++)
                    for (uint32_t a=0;a<n;a++)
                        for (uint32_t b=0;b<n;b++)
                            for (uint32_t c=0;c<n;c++)
                                temp_element1[i*n*n+b*n+c]+=element[a*n*n+b*n+c]*neg_L[a*n+i];
                for (uint32_t i=0;i<n;i++)
                    for (uint32_t b=0;b<n;b++)
                        for (uint32_t j=0;j<n;j++)
                            for (uint32_t c=0;c<n;c++)
                                temp_element2[i*n*n+j*n+c]+=temp_element1[i*n*n+b*n+c]*neg_L[b*n+j];
                memset(temp_element1.data(),0,sizeof(float)*n*n*n);
                for (uint32_t i=0;i<n;i++)
                    for (uint32_t j=0;j<n;j++)
                        for (uint32_t k=0;k<n;k++)
                            for (uint32_t c=0;c<n;c++)
                                temp_element1[i*n*n+j*n+k]+=temp_element2[i*n*n+j*n+c]*neg_L[c*n+k];
                quantize_and_overwrite2(temp_element1.data(),element_f_qb,n*n*n,param_eb,neg_2_param_eb,zero_freq_bound);
                memset(temp_element2.data(),0,sizeof(float)*n*n*n);
                for (uint32_t i=0;i<n;i++)
                    for (uint32_t j=0;j<n;j++)
                        for (uint32_t k=0;k<n;k++)
                            for (uint32_t c=0;c<n;c++)
                                temp_element2[i*n*n+j*n+c]+=temp_element1[i*n*n+j*n+k]*L[k*n+c];
                memset(temp_element1.data(),0,sizeof(float)*n*n*n);
                for (uint32_t i=0;i<n;i++)
                    for (uint32_t j=0;j<n;j++)
                        for (uint32_t b=0;b<n;b++)
                            for (uint32_t c=0;c<n;c++)
                                temp_element1[i*n*n+b*n+c]+=temp_element2[i*n*n+j*n+c]*L[j*n+b];
                memset(temp_element2.data(),0,sizeof(float)*n*n*n);
                for (uint32_t i=0;i<n;i++)
                    for (uint32_t a=0;a<n;a++)
                        for (uint32_t b=0;b<n;b++)
                            for (uint32_t c=0;c<n;c++)
                                temp_element2[a*n*n+b*n+c]+=temp_element1[i*n*n+b*n+c]*L[i*n+a];
                quantize_err(temp_element2.data(),element,element_err_qb,n*n*n,eb,neg_2_eb);
            }
            // 把三种qb无损压缩到一起
            vector<unsigned char> huffman_bs_vec((2*header.element_size+8)*(end_index-begin_index)*sizeof(int64_t));
            printf("begin_index= %lu, end_index= %lu\n",begin_index,end_index);
            unsigned char *huffman_bs=huffman_bs_vec.data();
            huffman_encoder_class<int64_t> corner_encoder;
            corner_encoder.preprocess_encode(corner_qb+8*begin_index,8*(end_index-begin_index));
            corner_encoder.save(huffman_bs);
            printf("huffman_bs_len= %lu\n",huffman_bs-huffman_bs_vec.data());
            corner_encoder.encode(corner_qb+8*begin_index,8*(end_index-begin_index),huffman_bs);
            printf("huffman_bs_len= %lu\n",huffman_bs-huffman_bs_vec.data());
            huffman_encoder_class<int64_t> f_encoder;
            f_encoder.preprocess_encode(f_qb+begin_index*header.element_size,header.element_size*(end_index-begin_index));
            f_encoder.save(huffman_bs);
            printf("huffman_bs_len= %lu\n",huffman_bs-huffman_bs_vec.data());
            f_encoder.encode(f_qb+begin_index*header.element_size,header.element_size*(end_index-begin_index),huffman_bs);
            printf("huffman_bs_len= %lu\n",huffman_bs-huffman_bs_vec.data());
            huffman_encoder_class<int64_t> err_encoder;
            err_encoder.preprocess_encode(err_qb+begin_index*header.element_size,header.element_size*(end_index-begin_index));
            err_encoder.save(huffman_bs);
            printf("huffman_bs_len= %lu\n",huffman_bs-huffman_bs_vec.data());
            err_encoder.encode(err_qb+begin_index*header.element_size,header.element_size*(end_index-begin_index),huffman_bs);
            printf("huffman_bs_len= %lu\n",huffman_bs-huffman_bs_vec.data());
            uint64_t huffman_bs_len=huffman_bs-huffman_bs_vec.data();
            uint64_t max_zstd_bs_len=ZSTD_compressBound(huffman_bs_len);
            zstd_bs[block_id].resize(max_zstd_bs_len);
            uint64_t zstd_bs_len=ZSTD_compress(zstd_bs[block_id].data(),max_zstd_bs_len,huffman_bs_vec.data(),huffman_bs_len,3);
            zstd_bs[block_id].resize(zstd_bs_len);
        }
        double end_time=omp_get_wtime();
        printf("Compression time= %f seconds\n",end_time-start_time);
        printf("Bandwidth per cpu= %f MB/s\n",element_num_to_process*header.element_size*header.word_size/((end_time-start_time)*1024*1024*available_cpus));
        uint64_t total_compressed_size=0;
        for (auto block_id:blocks_to_process)
            total_compressed_size+=zstd_bs[block_id].size();
        printf("Total compressed size= %lu bytes\n",total_compressed_size);
        float br=static_cast<float>(total_compressed_size*8)/(element_num_to_process*header.element_size);
        float cr=static_cast<float>(header.word_size*8)/br;
        printf("BR= %f bits/value, CR= %f\n",br,cr);
    }
}

#endif