#pragma once

namespace mrcz
{
    template<typename T>
    void param_class<T>::compress_mrcz(string field)
    {
        const bool experimental=method=="mrcz-experimental";
        uint32_t n0=header.element_shape[0];
        uint32_t n1=header.element_shape[1];
        uint32_t n2=header.element_shape[2];
        uint32_t n12=n1*n2;
        uint32_t n012=n0*n1*n2;
        vector<T> L0=L[0];
        vector<T> L1=L[1];
        vector<T> L2=L[2];
        vector<T> neg_L0=neg_L[0];
        vector<T> neg_L1=neg_L[1];
        vector<T> neg_L2=neg_L[2];
        T abs_eb=this->abs_eb[field];
        T coeff_eb=abs_eb*coeff_eb_scale;
        T zero_coeff_bound=coeff_eb*zero_coeff_qb_bound;
        T neg_2_coeff_eb=1/(2*coeff_eb);
        T neg_2_eb=1/(2*abs_eb);
        int64_t *coeff_qb=new int64_t[element_num_to_process*n012];
        int64_t *err_qb=new int64_t[element_num_to_process*n012];
        double start_time=omp_get_wtime();
        #pragma omp parallel for schedule(static) num_threads(max_cpus)
        for (uint64_t block_id=0;block_id<partitions.size();block_id++)
        {
            auto &block=partitions[block_id];
            uint64_t begin_index=block.begin_index;
            uint64_t end_index=block.end_index;
            // 然后量化剩余的数据
            vector<T> temp_element1(n012,0),temp_element2(n012,0);
            for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
            {
                T *element=data.data()+element_id*n012;
                int64_t *element_f_qb=coeff_qb+element_id*n012;
                int64_t *element_err_qb=err_qb+element_id*n012;
                // 正变换
                forward_transform(temp_element1,temp_element2,element,n0,n1,n2,n12,n012,L0,L1,L2);
                // 量化
                if (experimental and if_zero_freq_expanding)
                    quantize_and_overwrite_array2(temp_element1.data(),element_f_qb,n012,coeff_eb,neg_2_coeff_eb,zero_coeff_bound);
                else
                    quantize_and_overwrite_array1(temp_element1.data(),element_f_qb,n012,coeff_eb,neg_2_coeff_eb);
                // 反变换
                inverse_transform(temp_element1,temp_element2,n0,n1,n2,n12,n012,neg_L0,neg_L1,neg_L2);
                quantize_err_array(temp_element2.data(),element,element_err_qb,n012,abs_eb,neg_2_eb);
            }
            // 无损压缩
            vector<unsigned char> huffman_bs_vec(2*n012*(end_index-begin_index)*sizeof(int64_t),0);
            unsigned char *huffman_bs=huffman_bs_vec.data();
            if (experimental and if_zero_map)
                lossless_compress_coeff_qb_with_zero_map(begin_index,end_index,coeff_qb,huffman_bs,n012);
            else
                lossless_compress_coeff_qb(begin_index,end_index,coeff_qb,huffman_bs,n012);
            huffman_encoder_class<int64_t> err_encoder;
            err_encoder.preprocess_encode(err_qb+begin_index*n012,n012*(end_index-begin_index));
            err_encoder.save(huffman_bs);
            err_encoder.encode(err_qb+begin_index*n012,n012*(end_index-begin_index),huffman_bs);
            uint64_t huffman_bs_len=huffman_bs-huffman_bs_vec.data();
            uint64_t max_zstd_bs_len=ZSTD_compressBound(huffman_bs_len);
            block.zstd_bs.resize(max_zstd_bs_len);
            uint64_t zstd_bs_len=ZSTD_compress(block.zstd_bs.data(),max_zstd_bs_len,huffman_bs_vec.data(),huffman_bs_len,3);
            block.zstd_bs.resize(zstd_bs_len);

        }
        double end_time=omp_get_wtime();
        delete[] coeff_qb;
        delete[] err_qb;
        printf("Compression time= %f seconds\n",end_time-start_time);
        printf("\033[91mField= %s, compression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),element_num_to_process*n012*header.word_size/((end_time-start_time)*1024*1024*max_cpus));
    }

    template<typename T>
    inline void param_class<T>::forward_transform(vector<T> &temp_element1,vector<T> &temp_element2,T *element,uint32_t n0,uint32_t n1,uint32_t n2,uint32_t n12,uint32_t n012,
        vector<T> &L0,vector<T> &L1,vector<T> &L2)
    {
        memset(temp_element1.data(),0,sizeof(T)*n012);
        for (uint32_t a=0;a<n0;a++)
            for (uint32_t i=0;i<n0;i++)
                for (uint32_t bn2pc=0;bn2pc<n1*n2;bn2pc++)
                        temp_element1[i*n12+bn2pc]+=element[a*n12+bn2pc]*L0[a*n0+i];
        memset(temp_element2.data(),0,sizeof(T)*n012);
        for (uint32_t in12=0;in12<n012;in12+=n12)
            for (uint32_t b=0;b<n1;b++)
                for (uint32_t j=0;j<n1;j++)
                    for (uint32_t c=0;c<n2;c++)
                        temp_element2[in12+j*n2+c]+=temp_element1[in12+b*n2+c]*L1[b*n1+j];
        memset(temp_element1.data(),0,sizeof(T)*n012);
        for (uint32_t in12pjn2=0;in12pjn2<n012;in12pjn2+=n2)
            for (uint32_t c=0;c<n2;c++)
                for (uint32_t k=0;k<n2;k++)
                    temp_element1[in12pjn2+k]+=temp_element2[in12pjn2+c]*L2[c*n2+k];
    }

    template<typename T>
    inline void param_class<T>::inverse_transform(vector<T> &temp_element1,vector<T> &temp_element2,uint32_t n0,uint32_t n1,uint32_t n2,uint32_t n12,uint32_t n012,
        vector<T> &neg_L0,vector<T> &neg_L1,vector<T> &neg_L2)
    {
        memset(temp_element2.data(),0,sizeof(T)*n012);
        for (uint32_t in12pjn2=0;in12pjn2<n012;in12pjn2+=n2)
            for (uint32_t k=0;k<n2;k++)
                for (uint32_t c=0;c<n2;c++)
                    temp_element2[in12pjn2+c]+=temp_element1[in12pjn2+k]*neg_L2[k*n2+c];
        memset(temp_element1.data(),0,sizeof(T)*n012);
        for (uint32_t in12=0;in12<n012;in12+=n12)
            for (uint32_t j=0;j<n1;j++)
                for (uint32_t b=0;b<n1;b++)
                    for (uint32_t c=0;c<n2;c++)
                        temp_element1[in12+b*n2+c]+=temp_element2[in12+j*n2+c]*neg_L1[j*n1+b];
        memset(temp_element2.data(),0,sizeof(T)*n012);
        for (uint32_t i=0;i<n0;i++)
            for (uint32_t a=0;a<n0;a++)
                for (uint32_t bn2pc=0;bn2pc<n12;bn2pc++)
                    temp_element2[a*n12+bn2pc]+=temp_element1[i*n12+bn2pc]*neg_L0[i*n0+a];
    }
}
