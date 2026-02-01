#pragma once

#include "jzzz_param.hpp"

namespace jzzz
{
    void lossless_decompress_f_qb_with_zero_map(header_class &header,uint64_t begin_index,uint64_t end_index,
        int64_t *f_qb,int64_t *compact_f_qb,const unsigned char *&huffman_bs);

    template<typename T>
    void param_class<T>::decompress(string field)
    {
        uint32_t n=header.element_shape[0];
        uint32_t n2=n*n;
        uint32_t n3=n2*n;
        T abs_eb=this->abs_eb[field];
        T coeff_eb=abs_eb*coeff_eb_scale;
        int64_t *corner_qb=new int64_t[8*element_num_to_process];
        int64_t *f_qb=new int64_t[element_num_to_process*header.element_size];
        int64_t *compact_f_qb=new int64_t[element_num_to_process*header.element_size];
        int64_t *no_corner_err_qb=new int64_t[element_num_to_process*(header.element_size-8)];
        decompressed_data.resize(element_num_to_process*header.element_size);
        double start_time=omp_get_wtime();
        #pragma omp parallel for schedule(static) num_threads(max_cpus)
        for (uint64_t block_id=0;block_id<blocks.size();block_id++)
        {
            uint64_t begin_index=blocks[block_id].begin_index;
            uint64_t end_index=blocks[block_id].end_index;
            // 先无损解压
            uint64_t huffman_bs_len=ZSTD_getFrameContentSize(blocks[block_id].zstd_bs.data(),blocks[block_id].zstd_bs.size());
            vector<unsigned char> huffman_bs_vec(huffman_bs_len);
            if (ZSTD_decompress(huffman_bs_vec.data(),huffman_bs_len,blocks[block_id].zstd_bs.data(),blocks[block_id].zstd_bs.size())!=huffman_bs_len)
                throw runtime_error("ZSTD_decompress error in decompress()");
            const unsigned char *huffman_bs=huffman_bs_vec.data();
            huffman_encoder_class<int64_t> corner_encoder;
            corner_encoder.load(huffman_bs);
            corner_encoder.decode(huffman_bs,corner_qb+8*begin_index,8*(end_index-begin_index));
            lossless_decompress_f_qb_with_zero_map(header,begin_index,end_index,f_qb,compact_f_qb,huffman_bs);
            huffman_encoder_class<int64_t> err_encoder;
            err_encoder.load(huffman_bs);
            err_encoder.decode(huffman_bs,no_corner_err_qb+(header.element_size-8)*begin_index,(header.element_size-8)*(end_index-begin_index));
            
            // 重构每个element角上的数据
            for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
            {
                T *element=decompressed_data.data()+element_id*header.element_size;
                int64_t *element_corner_qb=corner_qb+8*element_id;
                element[0]=element_corner_qb[0]*2*abs_eb;
                element[n-1]=element_corner_qb[1]*2*abs_eb;
                element[(n-1)*n]=element_corner_qb[2]*2*abs_eb;
                element[n2-1]=element_corner_qb[3]*2*abs_eb;
                element[(n-1)*n2]=element_corner_qb[4]*2*abs_eb;
                element[(n-1)*(n2+1)]=element_corner_qb[5]*2*abs_eb;
                element[(n-1)*n*(n+1)]=element_corner_qb[6]*2*abs_eb;
                element[n3-1]=element_corner_qb[7]*2*abs_eb;
            }
            // 然后重构剩余的数据
            for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
            {
                T *element=decompressed_data.data()+element_id*header.element_size;
                int64_t *element_f_qb=f_qb+element_id*header.element_size;
                int64_t *element_no_corner_err_qb=no_corner_err_qb+element_id*(header.element_size-8);
                if (is_same<T,float>::value)
                {
                    float *float_element_base_weight=reinterpret_cast<float*>(element_base_weight.data());
                    __m256 vec_corner_numbers=_mm256_set_ps(element[n3-1],element[(n-1)*n*(n+1)],element[(n-1)*(n2+1)],element[(n-1)*n2],
                                                            element[n2-1],element[(n-1)*n],element[n-1],element[0]);
                    for (auto pos:no_corner_pos)
                    {
                        uint32_t bias=8*pos;
                        __m256 vec_element_base_weight=_mm256_loadu_ps(float_element_base_weight+bias);
                        __m256 vec_mult=_mm256_mul_ps(vec_corner_numbers,vec_element_base_weight);
                        __m256 sum1=_mm256_hadd_ps(vec_mult,vec_mult);
                        __m256 sum2=_mm256_hadd_ps(sum1,sum1);
                        __m128 low=_mm256_castps256_ps128(sum2);
                        __m128 high=_mm256_extractf128_ps(sum2,1);
                        __m128 sum3=_mm_add_ps(low,high);
                        float total_weight=_mm_cvtss_f32(sum3);
                        element[pos]+=total_weight;
                    }
                }
                else
                {
                    double *double_element_base_weight=reinterpret_cast<double*>(element_base_weight.data());
                    __m256d vec_corner_numbers_1=_mm256_set_pd(element[n3-1],element[n3-n],element[n3-n2+n-1],element[n3-n2]);
                    __m256d vec_corner_numbers_2=_mm256_set_pd(element[n2-1],element[n2-n],element[n-1],element[0]);
                    for (auto pos:no_corner_pos)
                    {
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
                        element[pos]+=total_weight;
                    }
                }
                vector<T> temp_element1(n3,0),temp_element2(n3,0);
                for (uint32_t i=0;i<n3;i++)
                    temp_element1[i]=element_f_qb[i]*2*coeff_eb;
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
                
                for (uint32_t j=1;j<n-1;j++)
                    element[j]+=(temp_element2[j]-element_no_corner_err_qb[j-1]*2*abs_eb);
                for (uint32_t j=n;j<n2-n;j++)
                    element[j]+=(temp_element2[j]-element_no_corner_err_qb[j-2]*2*abs_eb);
                for (uint32_t j=n2-n+1;j<n2-1;j++)
                    element[j]+=(temp_element2[j]-element_no_corner_err_qb[j-3]*2*abs_eb);
                for (uint32_t j=n2;j<n3-n2;j++)
                    element[j]+=(temp_element2[j]-element_no_corner_err_qb[j-4]*2*abs_eb);
                for (uint32_t j=n3-n2+1;j<n3-n2+n-1;j++)
                    element[j]+=(temp_element2[j]-element_no_corner_err_qb[j-5]*2*abs_eb);
                for (uint32_t j=n3-n2+n;j<n3-n;j++)
                    element[j]+=(temp_element2[j]-element_no_corner_err_qb[j-6]*2*abs_eb);
                for (uint32_t j=n3-n+1;j<n3-1;j++)
                    element[j]+=(temp_element2[j]-element_no_corner_err_qb[j-7]*2*abs_eb);
            }
        }
        double end_time=omp_get_wtime();
        delete[] corner_qb;
        delete[] f_qb;
        delete[] compact_f_qb;
        delete[] no_corner_err_qb;
        printf("Decompression time= %f seconds\n",end_time-start_time);
        printf("\033[91mBandwidth per cpu= %f MB/s\n\033[0m",element_num_to_process*header.element_size*header.word_size/((end_time-start_time)*1024*1024*max_cpus));
    }

    __attribute__((always_inline)) void lossless_decompress_f_qb_with_zero_map(header_class &header,uint64_t begin_index,uint64_t end_index,int64_t *f_qb,
        int64_t *compact_f_qb,const unsigned char *&huffman_bs)
    {
        bool f_qb_all_zero=false;
        uint32_t all_zero_num=0;
        vector<bool> f_qb_zero_map(header.element_size,false);
        vector<uint32_t> f_qb_compact_pos(header.element_size,0);
        if ((huffman_bs[0]&0xC0)==0)
        {
            f_qb_all_zero=true;
            f_qb_zero_map.assign(header.element_size,true);
            huffman_bs+=1;
        }
        else if ((huffman_bs[0]&0xC0)==0x80)
            huffman_bs+=1;
        else
        {
            uint32_t used_bits=2;
            for (uint32_t pos=0;pos<header.element_size;pos++)
            {
                if ((huffman_bs[0]>>(7-used_bits))&0x01)
                    f_qb_zero_map[pos]=true;
                used_bits++;
                if (used_bits==8)
                {
                    huffman_bs++;
                    used_bits=0;
                }
            }
            if (used_bits!=0)
                huffman_bs++;
        }
        for (uint32_t pos=0;pos<header.element_size;pos++)
        {
            if (f_qb_zero_map[pos]==true)
                all_zero_num++;
            if (pos==0)
                f_qb_compact_pos[pos]=0;
            else
                f_qb_compact_pos[pos]=f_qb_compact_pos[pos-1]+(f_qb_zero_map[pos]?0:1);
        }
        if (!f_qb_all_zero)
        {
            huffman_encoder_class<int64_t> f_encoder;
            f_encoder.load(huffman_bs);
            f_encoder.decode(huffman_bs,compact_f_qb+header.element_size*begin_index,(header.element_size-all_zero_num)*(end_index-begin_index));
            for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
            {
                auto compact_element_f_qb=compact_f_qb+begin_index*header.element_size+(element_id-begin_index)*(header.element_size-all_zero_num);
                auto original_element_f_qb=f_qb+element_id*header.element_size;
                for (uint32_t pos=0;pos<header.element_size;pos++)
                    if (f_qb_zero_map[pos]==false)
                        original_element_f_qb[pos]=compact_element_f_qb[f_qb_compact_pos[pos]];
                    else
                        original_element_f_qb[pos]=0;
            }
        }
        else
            memset(f_qb+begin_index*header.element_size,0,sizeof(int64_t)*header.element_size*(end_index-begin_index));
    }
}