#pragma once

namespace mrcz
{
    template<typename T>
    void param_class<T>::decompress_irregular(string field)
    {
        uint32_t n0=header.element_shape[0];
        uint32_t n1=header.element_shape[1];
        uint32_t n2=header.element_shape[2];
        uint32_t n12=n1*n2;
        uint32_t n012=n0*n1*n2;
        vector<T> L0=L[0];
        vector<T> L1=L[1];
        vector<T> L2=L[2];
        T abs_eb=this->abs_eb[field];
        T coeff_eb=abs_eb*coeff_eb_scale;
        int64_t *coeff_qb=new int64_t[element_num_to_process*n012];
        int64_t *err_qb=new int64_t[element_num_to_process*n012];
        decompressed_data.resize(element_num_to_process*n012);
        double start_time=omp_get_wtime();
        #pragma omp parallel for schedule(static) num_threads(max_cpus)
        for (uint64_t block_id=0;block_id<blocks.size();block_id++)
        {
            uint64_t begin_index=blocks[block_id].begin_index;
            uint64_t end_index=blocks[block_id].end_index;
            // 无损解压
            uint64_t huffman_bs_len=ZSTD_getFrameContentSize(blocks[block_id].zstd_bs.data(),blocks[block_id].zstd_bs.size());
            vector<unsigned char> huffman_bs_vec(huffman_bs_len);
            if (ZSTD_decompress(huffman_bs_vec.data(),huffman_bs_len,blocks[block_id].zstd_bs.data(),blocks[block_id].zstd_bs.size())!=huffman_bs_len)
                throw runtime_error("ZSTD_decompress error in decompress()");
            const unsigned char *huffman_bs=huffman_bs_vec.data();
            if (if_zero_map)
                lossless_decompress_coeff_qb_with_zero_map(begin_index,end_index,coeff_qb,huffman_bs,n012);
            else
                lossless_decompress_coeff_qb(begin_index,end_index,coeff_qb,huffman_bs,n012);
            huffman_encoder_class<int64_t> err_encoder;
            err_encoder.load(huffman_bs);
            err_encoder.decode(huffman_bs,err_qb+n012*begin_index,n012*(end_index-begin_index));
            
            // 然后重构剩余的数据
            vector<T> temp_element1(n012,0),temp_element2(n012,0);
            for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
            {
                T *element=decompressed_data.data()+element_id*n012;
                int64_t *element_coeff_qb=coeff_qb+element_id*n012;
                int64_t *element_err_qb=err_qb+element_id*n012;
                // 重构coeff
                for (uint32_t i=0;i<n012;i++)
                    temp_element1[i]=element_coeff_qb[i]*2*coeff_eb;
                // 反变换
                inverse_transform_irregular(temp_element1,temp_element2,n0,n1,n2,n12,n012,L0,L1,L2);
                for (uint32_t j=0;j<n012;j++)
                    element[j]+=(temp_element2[j]-element_err_qb[j]*2*abs_eb);
            }
        }
        double end_time=omp_get_wtime();
        delete[] err_qb;
        delete[] coeff_qb;
        printf("Decompression time= %f seconds\n",end_time-start_time);
        printf("\033[91mBandwidth per cpu= %f MB/s\033[0m\n",element_num_to_process*n012*header.word_size/((end_time-start_time)*1024*1024*max_cpus));
    }

    // template<typename T>
    // inline void param_class<T>::reconstruct_corner_number(int64_t *corner_qb,uint64_t begin_index,uint64_t end_index,uint32_t n2,uint32_t n12,uint32_t n012,T abs_eb)
    // {
    //     for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
    //     {
    //         T *element=decompressed_data.data()+element_id*n012;
    //         int64_t *element_corner_qb=corner_qb+8*element_id;
    //         element[0]=element_corner_qb[0]*2*abs_eb;
    //         element[n2-1]=element_corner_qb[1]*2*abs_eb;
    //         element[n12-n2]=element_corner_qb[2]*2*abs_eb;
    //         element[n12-1]=element_corner_qb[3]*2*abs_eb;
    //         element[n012-n12]=element_corner_qb[4]*2*abs_eb;
    //         element[n012-n12+n2-1]=element_corner_qb[5]*2*abs_eb;
    //         element[n012-n2]=element_corner_qb[6]*2*abs_eb;
    //         element[n012-1]=element_corner_qb[7]*2*abs_eb;
    //     }
    // }

    // template<>
    // inline void param_class<float>::add_element_base(float *element,uint32_t n2,uint32_t n12,uint32_t n012)
    // {
    //     __m256 vec_corner_numbers=_mm256_set_ps(element[n012-1],element[n012-n2],element[n012-n12+n2-1],element[n012-n12],
    //                                             element[n12-1],element[n12-n2],element[n2-1],element[0]);
    //     for (auto pos:no_corner_pos)
    //     {
    //         uint32_t bias=8*pos;
    //         __m256 vec_element_base_weight=_mm256_loadu_ps(element_base_weight.data()+bias);
    //         __m256 vec_mult=_mm256_mul_ps(vec_corner_numbers,vec_element_base_weight);
    //         __m256 sum1=_mm256_hadd_ps(vec_mult,vec_mult);
    //         __m256 sum2=_mm256_hadd_ps(sum1,sum1);
    //         __m128 low=_mm256_castps256_ps128(sum2);
    //         __m128 high=_mm256_extractf128_ps(sum2,1);
    //         __m128 sum3=_mm_add_ps(low,high);
    //         float total_weight=_mm_cvtss_f32(sum3);
    //         element[pos]+=total_weight;
    //     }
    // }

    // template<>
    // inline void param_class<double>::add_element_base(double *element,uint32_t n2,uint32_t n12,uint32_t n012)
    // {
    //     __m256d vec_corner_numbers_1=_mm256_set_pd(element[n012-1],element[n012-n2],element[n012-n12+n2-1],element[n012-n12]);
    //     __m256d vec_corner_numbers_2=_mm256_set_pd(element[n12-1],element[n12-n2],element[n2-1],element[0]);
    //     for (auto pos:no_corner_pos)
    //     {
    //         uint32_t bias=8*pos;
    //         __m256d vec_element_base_weight_1=_mm256_loadu_pd(element_base_weight.data()+bias+4);
    //         __m256d vec_element_base_weight_2=_mm256_loadu_pd(element_base_weight.data()+bias);
    //         __m256d vec_mult_1=_mm256_mul_pd(vec_corner_numbers_1,vec_element_base_weight_1);
    //         __m256d vec_mult_2=_mm256_mul_pd(vec_corner_numbers_2,vec_element_base_weight_2);
    //         __m256d vec_sum=_mm256_add_pd(vec_mult_1,vec_mult_2);
    //         __m128d low=_mm256_castpd256_pd128(vec_sum);
    //         __m128d high=_mm256_extractf128_pd(vec_sum,1);
    //         __m128d sum1=_mm_add_pd(low,high);
    //         __m128d sum2=_mm_hadd_pd(sum1,sum1);
    //         double total_weight=_mm_cvtsd_f64(sum2);
    //         element[pos]+=total_weight;
    //     }
    // }
}