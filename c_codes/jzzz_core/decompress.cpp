#ifndef __DECOMPRESS_CPP__
#define __DECOMPRESS_CPP__

#include <omp.h>
#include <immintrin.h>
#include "jzzz_param.hpp"
#include "huffman.cpp"
#include <cmath>
#include "../third_party/zstd/lib/zstd.h"

namespace jzzz
{
    void param_class::decompress()
    {
        const char *slurm_cpus=getenv("SLURM_CPUS_PER_TASK");
        int available_cpus=0;
        if (slurm_cpus!=nullptr)
            available_cpus=atoi(slurm_cpus);
        else
            throw runtime_error("没有获取到每个task可用的cpu数量");
        printf("Available cpus per task= %d\n",available_cpus);
        uint32_t n=header.element_shape[0];
        uint32_t n2=n*n;
        uint32_t n3=n2*n;
        param_eb=eb*param_eb_scale;
        int64_t *corner_qb=new int64_t[8*element_num_to_process];
        int64_t *f_qb=new int64_t[element_num_to_process*header.element_size];
        int64_t *no_corner_err_qb=new int64_t[element_num_to_process*(header.element_size-8)];
        decompressed_data.resize(element_num_to_process*header.element_size);
        double start_time=omp_get_wtime();
        #pragma omp parallel for schedule(static) num_threads(available_cpus)
        for (uint64_t fake_block_id=0;fake_block_id<blocks_to_process.size();fake_block_id++) //这个block_id是假的，指的是在需要处理的blocks中的编号，而不是header文件中的编号
        {
            uint64_t block_id=blocks_to_process[fake_block_id];
            // printf("CPU %d decompressing block %ld\n",omp_get_thread_num(),block_id);
            uint64_t begin_index=blocks[block_id].begin_index;
            uint64_t end_index=blocks[block_id].end_index;
            // 先无损解压
            uint64_t huffman_bs_len=ZSTD_getFrameContentSize(zstd_bs[block_id].data(),zstd_bs[block_id].size());
            vector<unsigned char> huffman_bs_vec(huffman_bs_len);
            if (ZSTD_decompress(huffman_bs_vec.data(),huffman_bs_len,zstd_bs[block_id].data(),zstd_bs[block_id].size())!=huffman_bs_len)
                throw runtime_error("ZSTD_decompress error in decompress()");
            const unsigned char *huffman_bs=huffman_bs_vec.data();
            huffman_encoder_class<int64_t> corner_encoder;
            corner_encoder.load(huffman_bs);
            corner_encoder.decode(huffman_bs,corner_qb+8*begin_index,8*(end_index-begin_index));
            huffman_encoder_class<int64_t> f_encoder;
            f_encoder.load(huffman_bs);
            f_encoder.decode(huffman_bs,f_qb+header.element_size*begin_index,header.element_size*(end_index-begin_index));
            huffman_encoder_class<int64_t> err_encoder;
            err_encoder.load(huffman_bs);
            err_encoder.decode(huffman_bs,no_corner_err_qb+(header.element_size-8)*begin_index,(header.element_size-8)*(end_index-begin_index));
            // 重构每个element角上的数据
            for (uint64_t i=begin_index;i<end_index;i++)
            {
                float *element=decompressed_data.data()+i*header.element_size;
                int64_t *element_corner_qb=corner_qb+8*i;
                element[0]=element_corner_qb[0]*2*eb;
                element[n-1]=element_corner_qb[1]*2*eb;
                element[(n-1)*n]=element_corner_qb[2]*2*eb;
                element[n*n-1]=element_corner_qb[3]*2*eb;
                element[(n-1)*n*n]=element_corner_qb[4]*2*eb;
                element[(n-1)*(n*n+1)]=element_corner_qb[5]*2*eb;
                element[(n-1)*n*(n+1)]=element_corner_qb[6]*2*eb;
                element[n*n*n-1]=element_corner_qb[7]*2*eb;
            }
            // 然后重构剩余的数据
            for (uint64_t i=begin_index;i<end_index;i++)
            {
                float *element=decompressed_data.data()+i*header.element_size;
                int64_t *element_f_qb=f_qb+i*header.element_size;
                int64_t *element_no_corner_err_qb=no_corner_err_qb+i*(header.element_size-8);
                __m256 vec_corner_numbers=_mm256_set_ps(element[n*n*n-1],element[(n-1)*n*(n+1)],element[(n-1)*(n*n+1)],element[(n-1)*n*n],
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
                            element[pos]+=total_weight;
                        }
                element[0]/=2;
                element[n-1]/=2;
                element[(n-1)*n]/=2;
                element[n*n-1]/=2;
                element[(n-1)*n*n]/=2;
                element[(n-1)*(n*n+1)]/=2;
                element[(n-1)*n*(n+1)]/=2;
                element[n*n*n-1]/=2;
                vector<float> temp_element1(n*n*n,0),temp_element2(n*n*n,0);
                for (uint32_t i=0;i<n*n*n;i++)
                    temp_element1[i]=element_f_qb[i]*2*param_eb;
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
                for (uint32_t j=1;j<n-1;j++)
                    element[j]+=(temp_element2[j]-element_no_corner_err_qb[j-1]*2*eb);
                for (uint32_t j=n;j<n2-n;j++)
                    element[j]+=(temp_element2[j]-element_no_corner_err_qb[j-2]*2*eb);
                for (uint32_t j=n2-n+1;j<n2-1;j++)
                    element[j]+=(temp_element2[j]-element_no_corner_err_qb[j-3]*2*eb);
                for (uint32_t j=n2;j<n3-n2;j++)
                    element[j]+=(temp_element2[j]-element_no_corner_err_qb[j-4]*2*eb);
                for (uint32_t j=n3-n2+1;j<n3-n2+n-1;j++)
                    element[j]+=(temp_element2[j]-element_no_corner_err_qb[j-5]*2*eb);
                for (uint32_t j=n3-n2+n;j<n3-n;j++)
                    element[j]+=(temp_element2[j]-element_no_corner_err_qb[j-6]*2*eb);
                for (uint32_t j=n3-n+1;j<n3-1;j++)
                    element[j]+=(temp_element2[j]-element_no_corner_err_qb[j-7]*2*eb);
            }
        }
        double end_time=omp_get_wtime();
        delete[] corner_qb;
        delete[] f_qb;
        delete[] no_corner_err_qb;
        printf("Decompression time= %f seconds\n",end_time-start_time);
        printf("Bandwidth per cpu= %f MB/s\n",element_num_to_process*header.element_size*header.word_size/((end_time-start_time)*1024*1024*available_cpus));
    }
}

#endif