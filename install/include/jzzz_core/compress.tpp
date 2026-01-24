#pragma once

#include "jzzz_param.hpp"

namespace jzzz
{
    template<typename T>
    void param_class<T>::compress(string field)
    {
        uint32_t n=header.element_shape[0];
        uint32_t n2=n*n;
        uint32_t n3=n2*n;
        T abs_eb=this->abs_eb[field];
        T coeff_eb=abs_eb*coeff_eb_scale;
        T zero_freq_bound=coeff_eb*zero_f_qb_bound;
        T neg_2_coeff_eb=1/(2*coeff_eb);
        T neg_2_eb=1/(2*abs_eb);
        int64_t *corner_qb=new int64_t[8*element_num_to_process];
        int64_t *f_qb=new int64_t[element_num_to_process*header.element_size];
        int64_t *err_qb=new int64_t[element_num_to_process*header.element_size];
        int64_t *no_corner_err_qb=new int64_t[element_num_to_process*(header.element_size-8)];
        printf("Begin compressing\n");
        double start_time=omp_get_wtime();
        #pragma omp parallel for schedule(static) num_threads(max_cpus)
        for (uint64_t block_id=0;block_id<blocks.size();block_id++)
        {
            auto &block=blocks[block_id];
            uint64_t begin_index=block.begin_index;
            uint64_t end_index=block.end_index;
            // 先量化位于每个element位于角上的点，覆盖原数据
            // 之后可能要改成压缩率更高的方法
            vector<T> corner_numbers(8*(end_index-begin_index));
            for (uint64_t i=begin_index;i<end_index;i++)
            {
                T *element=data.data()+i*header.element_size;
                int64_t *element_corner_qb=corner_qb+8*i;
                quantize_and_overwrite(element[0],element_corner_qb[0],abs_eb,neg_2_eb);
                quantize_and_overwrite(element[n-1],element_corner_qb[1],abs_eb,neg_2_eb);
                quantize_and_overwrite(element[n2-n],element_corner_qb[2],abs_eb,neg_2_eb);
                quantize_and_overwrite(element[n2-1],element_corner_qb[3],abs_eb,neg_2_eb);
                quantize_and_overwrite(element[n3-n2],element_corner_qb[4],abs_eb,neg_2_eb);
                quantize_and_overwrite(element[n3-n2+n-1],element_corner_qb[5],abs_eb,neg_2_eb);
                quantize_and_overwrite(element[n3-n],element_corner_qb[6],abs_eb,neg_2_eb);
                quantize_and_overwrite(element[n3-1],element_corner_qb[7],abs_eb,neg_2_eb);
            }
            // 然后量化剩余的数据
            for (uint64_t i=begin_index;i<end_index;i++)
            {
                T *element=data.data()+i*header.element_size;
                int64_t *element_f_qb=f_qb+i*header.element_size;
                int64_t *element_err_qb=err_qb+i*header.element_size;
                // 使用8个角上的数据插值其他数据
                if (is_same<T,float>::value)
                {
                    float *float_element_base_weight=reinterpret_cast<float*>(element_base_weight.data());
                    __m256 vec_corner_numbers=_mm256_set_ps(element[n3-1],element[n3-n],element[n3-n2+n-1],element[n3-n2],
                                                            element[n2-1],element[n2-n],element[n-1],element[0]);
                    for (uint32_t i0=0;i0<n;i0++)
                        for (uint32_t i1=0;i1<n;i1++)
                            for (uint32_t i2=0;i2<n;i2++)
                            {
                                uint32_t pos=i0*n2+i1*n+i2;
                                uint32_t bias=8*pos;
                                __m256 vec_element_base_weight=_mm256_loadu_ps(float_element_base_weight+bias);
                                __m256 vec_mult=_mm256_mul_ps(vec_corner_numbers,vec_element_base_weight);
                                __m256 sum1=_mm256_hadd_ps(vec_mult,vec_mult);
                                __m256 sum2=_mm256_hadd_ps(sum1,sum1);
                                __m128 low=_mm256_castps256_ps128(sum2);
                                __m128 high=_mm256_extractf128_ps(sum2,1);
                                __m128 sum3=_mm_add_ps(low,high);
                                float total_weight=_mm_cvtss_f32(sum3);
                                element[pos]-=total_weight;
                            }
                }
                else
                {
                    double *double_element_base_weight=reinterpret_cast<double*>(element_base_weight.data());
                    __m256d vec_corner_numbers_1=_mm256_set_pd(element[n3-1],element[n3-n],element[n3-n2+n-1],element[n3-n2]);
                    __m256d vec_corner_numbers_2=_mm256_set_pd(element[n2-1],element[n2-n],element[n-1],element[0]);
                    for (uint32_t i0=0;i0<n;i0++)
                        for (uint32_t i1=0;i1<n;i1++)
                            for (uint32_t i2=0;i2<n;i2++)
                            {
                                uint32_t pos=i0*n2+i1*n+i2;
                                uint32_t bias=8*pos;
                                __m256d vec_element_base_weight_1=_mm256_loadu_pd(double_element_base_weight+bias+4);
                                __m256d vec_element_base_weight_2=_mm256_loadu_pd(double_element_base_weight+bias);
                                __m256d vec_mult_1=_mm256_mul_pd(vec_corner_numbers_1,vec_element_base_weight_1);
                                __m256d vec_mult_2=_mm256_mul_pd(vec_corner_numbers_2,vec_element_base_weight_2);
                                __m256d vec_sum=_mm256_add_pd(vec_mult_1,vec_mult_2);
                                __m128d low=_mm256_castpd256_pd128(vec_sum);
                                __m128d high=_mm256_extractf128_pd(vec_sum,1);
                                __m128d sum1=_mm_add_pd(low,high);
                                __m128d sum2=_mm_hadd_pd(sum1,sum1);
                                double total_weight=_mm_cvtsd_f64(sum2);
                                element[pos]-=total_weight;
                            }
                }
                // 正变换
                vector<T> temp_element1(n3,0),temp_element2(n3,0);
                for (uint32_t i=0;i<n;i++)
                    for (uint32_t a=0;a<n;a++)
                        for (uint32_t b=0;b<n;b++)
                            for (uint32_t c=0;c<n;c++)
                                temp_element1[i*n2+b*n+c]+=element[a*n2+b*n+c]*neg_L[a*n+i];
                for (uint32_t i=0;i<n;i++)
                    for (uint32_t b=0;b<n;b++)
                        for (uint32_t j=0;j<n;j++)
                            for (uint32_t c=0;c<n;c++)
                                temp_element2[i*n2+j*n+c]+=temp_element1[i*n2+b*n+c]*neg_L[b*n+j];
                memset(temp_element1.data(),0,sizeof(T)*n3);
                for (uint32_t i=0;i<n;i++)
                    for (uint32_t j=0;j<n;j++)
                        for (uint32_t k=0;k<n;k++)
                            for (uint32_t c=0;c<n;c++)
                                temp_element1[i*n2+j*n+k]+=temp_element2[i*n2+j*n+c]*neg_L[c*n+k];
                quantize_and_overwrite_array2(temp_element1.data(),element_f_qb,n3,coeff_eb,neg_2_coeff_eb,zero_freq_bound);
                // 反变换
                memset(temp_element2.data(),0,sizeof(T)*n3);
                for (uint32_t i=0;i<n;i++)
                    for (uint32_t j=0;j<n;j++)
                        for (uint32_t k=0;k<n;k++)
                            for (uint32_t c=0;c<n;c++)
                                temp_element2[i*n2+j*n+c]+=temp_element1[i*n2+j*n+k]*L[k*n+c];
                memset(temp_element1.data(),0,sizeof(T)*n3);
                for (uint32_t i=0;i<n;i++)
                    for (uint32_t j=0;j<n;j++)
                        for (uint32_t b=0;b<n;b++)
                            for (uint32_t c=0;c<n;c++)
                                temp_element1[i*n2+b*n+c]+=temp_element2[i*n2+j*n+c]*L[j*n+b];
                memset(temp_element2.data(),0,sizeof(T)*n3);
                for (uint32_t i=0;i<n;i++)
                    for (uint32_t a=0;a<n;a++)
                        for (uint32_t b=0;b<n;b++)
                            for (uint32_t c=0;c<n;c++)
                                temp_element2[a*n2+b*n+c]+=temp_element1[i*n2+b*n+c]*L[i*n+a];
                quantize_err_array(temp_element2.data(),element,element_err_qb,n3,abs_eb,neg_2_eb);
            }
            //去掉err_qb中位于角上的qb
            for (uint64_t i=begin_index;i<end_index;i++)
            {
                int64_t *element_err_qb=err_qb+i*header.element_size;
                int64_t *element_no_corner_err_qb=no_corner_err_qb+i*(header.element_size-8);
                memcpy(element_no_corner_err_qb,element_err_qb+1,sizeof(int64_t)*(n-2));
                memcpy(element_no_corner_err_qb+n-2,element_err_qb+n,sizeof(int64_t)*(n2-2*n));
                memcpy(element_no_corner_err_qb+n2-n-2,element_err_qb+n2-n+1,sizeof(int64_t)*(n-2));
                memcpy(element_no_corner_err_qb+n2-4,element_err_qb+n2,sizeof(int64_t)*(n3-2*n2));
                memcpy(element_no_corner_err_qb+n3-n2-4,element_err_qb+n3-n2+1,sizeof(int64_t)*(n-2));
                memcpy(element_no_corner_err_qb+n3-n2+n-6,element_err_qb+n3-n2+n,sizeof(int64_t)*(n2-2*n));
                memcpy(element_no_corner_err_qb+n3-n-6,element_err_qb+n3-n+1,sizeof(int64_t)*(n-2));
            }
            // 把三种qb无损压缩到一起
            vector<unsigned char> huffman_bs_vec((2*header.element_size+8)*(end_index-begin_index)*sizeof(int64_t));
            unsigned char *huffman_bs=huffman_bs_vec.data();
            huffman_encoder_class<int64_t> corner_encoder;
            corner_encoder.preprocess_encode(corner_qb+8*begin_index,8*(end_index-begin_index));
            corner_encoder.save(huffman_bs);
            corner_encoder.encode(corner_qb+8*begin_index,8*(end_index-begin_index),huffman_bs);
            huffman_encoder_class<int64_t> f_encoder;
            f_encoder.preprocess_encode(f_qb+begin_index*header.element_size,header.element_size*(end_index-begin_index));
            f_encoder.save(huffman_bs);
            f_encoder.encode(f_qb+begin_index*header.element_size,header.element_size*(end_index-begin_index),huffman_bs);
            huffman_encoder_class<int64_t> err_encoder;
            err_encoder.preprocess_encode(no_corner_err_qb+begin_index*(header.element_size-8),(header.element_size-8)*(end_index-begin_index));
            err_encoder.save(huffman_bs);
            err_encoder.encode(no_corner_err_qb+begin_index*(header.element_size-8),(header.element_size-8)*(end_index-begin_index),huffman_bs);
            uint64_t huffman_bs_len=huffman_bs-huffman_bs_vec.data();
            uint64_t max_zstd_bs_len=ZSTD_compressBound(huffman_bs_len);
            block.zstd_bs.resize(max_zstd_bs_len);
            uint64_t zstd_bs_len=ZSTD_compress(block.zstd_bs.data(),max_zstd_bs_len,huffman_bs_vec.data(),huffman_bs_len,3);
            block.zstd_bs.resize(zstd_bs_len);
        }
        double end_time=omp_get_wtime();
        delete[] corner_qb;
        delete[] f_qb;
        delete[] err_qb;
        delete[] no_corner_err_qb;
        printf("Compression time= %f seconds\n",end_time-start_time);
        printf("\033[91mBandwidth per cpu= %f MB/s\n\033[0m",element_num_to_process*header.element_size*header.word_size/((end_time-start_time)*1024*1024*max_cpus));
    }
}