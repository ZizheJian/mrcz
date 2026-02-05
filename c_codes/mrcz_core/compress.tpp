#pragma once

namespace mrcz
{
    template<typename T>
    void param_class<T>::compress(string field)
    {
        uint32_t n=header.element_shape[0];
        uint32_t n2=n*n;
        uint32_t n3=n2*n;
        T abs_eb=this->abs_eb[field];
        T coeff_eb=abs_eb*coeff_eb_scale;
        T zero_coeff_bound=coeff_eb*zero_coeff_qb_bound;
        T neg_2_coeff_eb=1/(2*coeff_eb);
        T neg_2_eb=1/(2*abs_eb);
        int64_t *corner_qb=new int64_t[8*element_num_to_process];
        int64_t *coeff_qb=new int64_t[element_num_to_process*n3];
        int64_t *err_qb=new int64_t[element_num_to_process*n3];
        int64_t *no_corner_err_qb=new int64_t[element_num_to_process*(n3-8)];
        // ofstream coeff_file(data_path+"_coeff",ios::binary);
        // ofstream err_file(data_path+"_err",ios::binary);
        printf("Begin compressing\n");
        double start_time=omp_get_wtime();
        #pragma omp parallel for schedule(static) num_threads(max_cpus)
        for (uint64_t block_id=0;block_id<blocks.size();block_id++)
        {
            auto &block=blocks[block_id];
            uint64_t begin_index=block.begin_index;
            uint64_t end_index=block.end_index;
            // 先量化位于每个element位于角上的点，覆盖原数据
            // 之后可能要改成压缩率更高的方法。因为可能使用跨element预测，所以放在element级for循环外面
            if (if_corner_based_prediction)
                quantize_and_overwrite_corner_numbers(corner_qb,begin_index,end_index,n,n2,n3,abs_eb,neg_2_eb);
            // 然后量化剩余的数据
            vector<T> temp_element1(n3,0),temp_element2(n3,0);
            for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
            {
                T *element=data.data()+element_id*n3;
                int64_t *element_f_qb=coeff_qb+element_id*n3;
                int64_t *element_err_qb=err_qb+element_id*n3;
                int64_t *element_no_corner_err_qb=no_corner_err_qb+element_id*(n3-8);
                // 使用8个角上的数据插值其他数据
                if (if_corner_based_prediction)
                    subtract_element_base(element,n,n2,n3);
                // 正变换
                forward_transform(temp_element1,temp_element2,element,n,n2,n3);
                // coeff_file.write(reinterpret_cast<char*>(temp_element1.data()),sizeof(T)*n3);
                // 量化
                if (if_zero_freq_expanding)
                    quantize_and_overwrite_array2(temp_element1.data(),element_f_qb,n3,coeff_eb,neg_2_coeff_eb,zero_coeff_bound);
                else
                    quantize_and_overwrite_array1(temp_element1.data(),element_f_qb,n3,coeff_eb,neg_2_coeff_eb);
                // 反变换
                inverse_transform(temp_element1,temp_element2,n,n2,n3);
                // vector<T> err(n3,0);
                // for (uint32_t i=0;i<n3;i++)
                //     err[i]=temp_element2[i]-element[i];
                // err_file.write(reinterpret_cast<char*>(err.data()),sizeof(T)*n3);
                quantize_err_array(temp_element2.data(),element,element_err_qb,n3,abs_eb,neg_2_eb);
                if (if_corner_based_prediction)
                    convert_err_eb_to_no_corner_err_eb(element_no_corner_err_qb,element_err_qb,n,n2,n3);
            }
            // 无损压缩
            vector<unsigned char> huffman_bs_vec((2*n3+8)*(end_index-begin_index)*sizeof(int64_t),0);
            unsigned char *huffman_bs=huffman_bs_vec.data();
            if (if_corner_based_prediction)
            {
                huffman_encoder_class<int64_t> corner_encoder;
                corner_encoder.preprocess_encode(corner_qb+8*begin_index,8*(end_index-begin_index));
                corner_encoder.save(huffman_bs);
                corner_encoder.encode(corner_qb+8*begin_index,8*(end_index-begin_index),huffman_bs);
            }
            if (if_zero_map)
                lossless_compress_coeff_qb_with_zero_map(begin_index,end_index,coeff_qb,huffman_bs);
            else
                lossless_compress_coeff_qb(begin_index,end_index,coeff_qb,huffman_bs);
            if (if_corner_based_prediction)
            {
                huffman_encoder_class<int64_t> err_encoder;
                err_encoder.preprocess_encode(no_corner_err_qb+begin_index*(n3-8),(n3-8)*(end_index-begin_index));
                err_encoder.save(huffman_bs);
                err_encoder.encode(no_corner_err_qb+begin_index*(n3-8),(n3-8)*(end_index-begin_index),huffman_bs);
            }
            else
            {
                huffman_encoder_class<int64_t> err_encoder;
                err_encoder.preprocess_encode(err_qb+begin_index*n3,n3*(end_index-begin_index));
                err_encoder.save(huffman_bs);
                err_encoder.encode(err_qb+begin_index*n3,n3*(end_index-begin_index),huffman_bs);
            }
            uint64_t huffman_bs_len=huffman_bs-huffman_bs_vec.data();
            uint64_t max_zstd_bs_len=ZSTD_compressBound(huffman_bs_len);
            block.zstd_bs.resize(max_zstd_bs_len);
            uint64_t zstd_bs_len=ZSTD_compress(block.zstd_bs.data(),max_zstd_bs_len,huffman_bs_vec.data(),huffman_bs_len,3);
            block.zstd_bs.resize(zstd_bs_len);
        }
        double end_time=omp_get_wtime();
        delete[] corner_qb;
        delete[] coeff_qb;
        delete[] err_qb;
        delete[] no_corner_err_qb;
        printf("Compression time= %f seconds\n",end_time-start_time);
        printf("\033[91mBandwidth per cpu= %f MB/s\n\033[0m",element_num_to_process*n3*header.word_size/((end_time-start_time)*1024*1024*max_cpus));
    }

    template<typename T>
    inline void param_class<T>::quantize_and_overwrite_corner_numbers(int64_t *corner_qb,uint64_t begin_index,uint64_t end_index,uint32_t n,uint32_t n2,uint32_t n3,
        T abs_eb,T neg_2_eb)
    {
        for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
        {
            T *element=data.data()+element_id*n3;
            int64_t *element_corner_qb=corner_qb+8*element_id;
            quantize_and_overwrite(element[0],element_corner_qb[0],abs_eb,neg_2_eb);
            quantize_and_overwrite(element[n-1],element_corner_qb[1],abs_eb,neg_2_eb);
            quantize_and_overwrite(element[n2-n],element_corner_qb[2],abs_eb,neg_2_eb);
            quantize_and_overwrite(element[n2-1],element_corner_qb[3],abs_eb,neg_2_eb);
            quantize_and_overwrite(element[n3-n2],element_corner_qb[4],abs_eb,neg_2_eb);
            quantize_and_overwrite(element[n3-n2+n-1],element_corner_qb[5],abs_eb,neg_2_eb);
            quantize_and_overwrite(element[n3-n],element_corner_qb[6],abs_eb,neg_2_eb);
            quantize_and_overwrite(element[n3-1],element_corner_qb[7],abs_eb,neg_2_eb);
        }
    }

    template<>
    inline void param_class<float>::subtract_element_base(float *element,uint32_t n,uint32_t n2,uint32_t n3)
    {
        __m256 vec_corner_numbers=_mm256_set_ps(element[n3-1],element[n3-n],element[n3-n2+n-1],element[n3-n2],
                                                element[n2-1],element[n2-n],element[n-1],element[0]);
        for (uint32_t pos=0;pos<n3;pos++)
        {
            uint32_t bias=8*pos;
            __m256 vec_element_base_weight=_mm256_loadu_ps(element_base_weight.data()+bias);
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

    template<>
    inline void param_class<double>::subtract_element_base(double *element,uint32_t n,uint32_t n2,uint32_t n3)
    {
        __m256d vec_corner_numbers_1=_mm256_set_pd(element[n3-1],element[n3-n],element[n3-n2+n-1],element[n3-n2]);
        __m256d vec_corner_numbers_2=_mm256_set_pd(element[n2-1],element[n2-n],element[n-1],element[0]);
        for (uint32_t pos=0;pos<n3;pos++)
        {
            uint32_t bias=8*pos;
            __m256d vec_element_base_weight_1=_mm256_loadu_pd(element_base_weight.data()+bias+4);
            __m256d vec_element_base_weight_2=_mm256_loadu_pd(element_base_weight.data()+bias);
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

    template<typename T>
    inline void param_class<T>::forward_transform(vector<T> &temp_element1,vector<T> &temp_element2,T *element,uint32_t n,uint32_t n2,uint32_t n3)
    {
        memset(temp_element1.data(),0,sizeof(T)*n3);
        for (uint32_t i=0;i<n;i++)
            for (uint32_t a=0;a<n;a++)
                for (uint32_t b=0;b<n;b++)
                    for (uint32_t c=0;c<n;c++)
                        temp_element1[i*n2+b*n+c]+=element[a*n2+b*n+c]*neg_L[a*n+i];
        memset(temp_element2.data(),0,sizeof(T)*n3);
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
    }

    template<typename T>
    inline void param_class<T>::inverse_transform(vector<T> &temp_element1,vector<T> &temp_element2,uint32_t n,uint32_t n2,uint32_t n3)
    {
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
    }

    template<typename T>
    inline void param_class<T>::convert_err_eb_to_no_corner_err_eb(int64_t *element_no_corner_err_qb,int64_t *element_err_qb,uint32_t n,uint32_t n2,uint32_t n3)
    {
        memcpy(element_no_corner_err_qb,element_err_qb+1,sizeof(int64_t)*(n-2));
        memcpy(element_no_corner_err_qb+n-2,element_err_qb+n,sizeof(int64_t)*(n2-2*n));
        memcpy(element_no_corner_err_qb+n2-n-2,element_err_qb+n2-n+1,sizeof(int64_t)*(n-2));
        memcpy(element_no_corner_err_qb+n2-4,element_err_qb+n2,sizeof(int64_t)*(n3-2*n2));
        memcpy(element_no_corner_err_qb+n3-n2-4,element_err_qb+n3-n2+1,sizeof(int64_t)*(n-2));
        memcpy(element_no_corner_err_qb+n3-n2+n-6,element_err_qb+n3-n2+n,sizeof(int64_t)*(n2-2*n));
        memcpy(element_no_corner_err_qb+n3-n-6,element_err_qb+n3-n+1,sizeof(int64_t)*(n-2));
    }

    template<typename T>
    inline void param_class<T>::lossless_compress_coeff_qb_with_zero_map(uint64_t begin_index,uint64_t end_index,int64_t *coeff_qb,unsigned char *&huffman_bs)
    {
        uint32_t n3=header.element_size;
        bool coeff_qb_all_zero=true;//coeff_qb全0
        bool coeff_qb_none_zero=true;//所有沿着element_id维度的柱子上至少有一个非0
        uint32_t all_zero_num=0;//沿着element_id的维度上全0的柱子数量
        vector<bool> coeff_qb_zero_map(n3,true);//沿着element_id的维度上全0的柱子的标记
        vector<uint32_t> coeff_qb_compact_pos(n3,0);//从非全0的柱子位置到紧凑存储时的位置的映射
        for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
            for (uint32_t pos=0;pos<n3;pos++)
                if (coeff_qb[element_id*n3+pos]!=0)
                    coeff_qb_zero_map[pos]=false;
        for (uint32_t pos=0;pos<n3;pos++)
        {
            coeff_qb_compact_pos[pos]=pos-all_zero_num;
            if (coeff_qb_zero_map[pos]==true)
            {
                coeff_qb_none_zero=false;
                all_zero_num++;
            }
            else
                coeff_qb_all_zero=false;
        }
        for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
        {
            auto original_element_f_qb=coeff_qb+element_id*n3;
            auto compact_element_f_qb=coeff_qb+begin_index*n3+(element_id-begin_index)*(n3-all_zero_num);
            for (uint32_t pos=0;pos<n3;pos++)
                if (coeff_qb_zero_map[pos]==false)
                    compact_element_f_qb[coeff_qb_compact_pos[pos]]=original_element_f_qb[pos];
        }
        
        if (coeff_qb_all_zero)
        {
            huffman_bs[0]=0;
            huffman_bs+=1;
        }
        else if (coeff_qb_none_zero)
        {
            huffman_bs[0]=0x80;
            huffman_bs+=1;
        }
        else
        {
            huffman_bs[0]=0x40;
            uint32_t used_bits=2;
            for (uint32_t pos=0;pos<n3;pos++)
            {
                if (coeff_qb_zero_map[pos]==true)
                    huffman_bs[0]|=((unsigned char)1)<<(7-used_bits);
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

        if (!coeff_qb_all_zero)
            if (if_coefficient_transpose)
            {
                int64_t *temp_coeff_qb=new int64_t[(n3-all_zero_num)*(end_index-begin_index)];
                for (uint32_t pos=0;pos<n3-all_zero_num;pos++)
                    for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
                        temp_coeff_qb[pos*(end_index-begin_index)+(element_id-begin_index)]=coeff_qb[begin_index*n3+(element_id-begin_index)*(n3-all_zero_num)+pos];
                huffman_encoder_class<int64_t> f_encoder;
                f_encoder.preprocess_encode(temp_coeff_qb,(n3-all_zero_num)*(end_index-begin_index));
                f_encoder.save(huffman_bs);
                f_encoder.encode(temp_coeff_qb,(n3-all_zero_num)*(end_index-begin_index),huffman_bs);
                delete[] temp_coeff_qb;
            }
            else
            {
                huffman_encoder_class<int64_t> f_encoder;
                f_encoder.preprocess_encode(coeff_qb+begin_index*n3,(n3-all_zero_num)*(end_index-begin_index));
                f_encoder.save(huffman_bs);
                f_encoder.encode(coeff_qb+begin_index*n3,(n3-all_zero_num)*(end_index-begin_index),huffman_bs);
            }
    }

    template<typename T>
    inline void param_class<T>::lossless_compress_coeff_qb(uint64_t begin_index,uint64_t end_index,int64_t *coeff_qb,unsigned char *&huffman_bs)
    {
        uint32_t n3=header.element_size;
        if (if_coefficient_transpose)
        {
            int64_t *temp_coeff_qb=new int64_t[n3*(end_index-begin_index)];
            for (uint32_t pos=0;pos<n3;pos++)
                for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
                    temp_coeff_qb[pos*(end_index-begin_index)+(element_id-begin_index)]=coeff_qb[element_id*n3+pos];
            huffman_encoder_class<int64_t> f_encoder;
            f_encoder.preprocess_encode(temp_coeff_qb,n3*(end_index-begin_index));
            f_encoder.save(huffman_bs);
            f_encoder.encode(temp_coeff_qb,n3*(end_index-begin_index),huffman_bs);
            delete[] temp_coeff_qb;
        }
        else
        {
            huffman_encoder_class<int64_t> f_encoder;
            f_encoder.preprocess_encode(coeff_qb+begin_index*n3,n3*(end_index-begin_index));
            f_encoder.save(huffman_bs);
            f_encoder.encode(coeff_qb+begin_index*n3,n3*(end_index-begin_index),huffman_bs);
        }
    }
}