#ifndef __COMPRESS_CPP__
#define __COMPRESS_CPP__

#include <omp.h>
#include <immintrin.h>
#include "jzzz_param.hpp"

namespace jzzz
{
    void param_class::compress()
    {
        uint32_t n=header.element_shape[0];
        param_eb=eb*param_eb_scale;
        float zero_freq_bound=param_eb*zero_f_qb_bound;
        float neg_2_param_eb=1/(2*param_eb);
        float neg_2_eb=1/(2*eb);
        f_qb.resize(element_num_to_process*header.element_size);
        qb.resize(element_num_to_process*header.element_size);
        double start_time=omp_get_wtime();
        #pragma omp parallel for schedule(static) num_threads(128)
        for (uint64_t block_id=0;block_id<blocks.size();block_id++) //这个block_id是假的，指的是在需要处理的blocks中的编号，而不是header文件中的编号
        {
            uint64_t begin_index=blocks[block_id].begin_index;
            uint64_t end_index=blocks[block_id].end_index;
            for (uint64_t i=begin_index;i<end_index;i++)
            {
                float *element=data.data()+i*header.element_size;
                int64_t *f_qb=this->f_qb.data()+i*header.element_size;
                int64_t *qb=this->qb.data()+i*header.element_size;
                __m256 vec_corner_numbers=_mm256_set_ps(element[n*n*n-1],element[(n-1)*n*n+(n-1)*n],element[(n-1)*n*n+n-1],element[(n-1)*n*n],
                                                    element[n*n-1],element[(n-1)*n],element[(n-1)*n*n],element[0]);
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
                for (uint32_t i=0;i<n*n*n;i++)
                {
                    if (-zero_freq_bound<=temp_element1[i] && temp_element1[i]<=zero_freq_bound)
                        f_qb[i]=0;
                    else
                        f_qb[i]=static_cast<int64_t>(temp_element1[i]*neg_2_param_eb);
                    temp_element1[i]=f_qb[i]*2*param_eb;
                }
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
                for (uint32_t i=0;i<n*n*n;i++)
                    qb[i]=static_cast<int64_t>((temp_element2[i]-element[i])*neg_2_eb);
            }
            
        }
        double end_time=omp_get_wtime();
        printf("Compression time= %f seconds\n",end_time-start_time);
        printf("Bandwidth per cpu= %f MB/s\n",element_num_to_process*header.element_size*header.word_size/((end_time-start_time)*1024*1024*128));
    }
}

#endif