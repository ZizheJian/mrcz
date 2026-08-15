#pragma once

namespace mrcz
{
    template<typename T>
    void param_class<T>::compress_regular(string field)
    {
        assert((header.element_shape[0]==header.element_shape[1]) and (header.element_shape[1]==header.element_shape[2]));
        uint32_t n0=header.element_shape[0];
        uint32_t n00=n0*n0;
        uint32_t n000=n00*n0;
        vector<T> L0=L[0];
        vector<T> neg_L0=neg_L[0];
        T abs_eb=this->abs_eb[field];
        T coeff_eb=abs_eb*coeff_eb_scale;
        T zero_coeff_bound=coeff_eb*zero_coeff_qb_bound;
        T neg_2_eb=1/(2*abs_eb);
        T neg_2_coeff_eb=1/(2*coeff_eb);
        int64_t *coeff_qb=new int64_t[element_num_to_process*n000];
        int64_t *err_qb=new int64_t[element_num_to_process*n000];
        // vector<double> coeff_br(max_cpus,0),err_br(max_cpus,0);
        double start_time=omp_get_wtime();
        #pragma omp parallel for schedule(static) num_threads(max_cpus)
        for (uint64_t block_id=0;block_id<partitions.size();block_id++)
        {
            auto &block=partitions[block_id];
            uint64_t begin_index=block.begin_index;
            uint64_t end_index=block.end_index;
            // 然后量化剩余的数据
            vector<T> temp_element1(n000,0),temp_element2(n000,0);
            vector<T> mean_abs_coeff(n000,0);
            for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
            {
                T *element=data.data()+element_id*n000;
                int64_t *element_f_qb=coeff_qb+element_id*n000;
                int64_t *element_err_qb=err_qb+element_id*n000;
                // 正变换
                forward_transform_regular(temp_element1,temp_element2,element,n0,n00,n000,L0);
                for (uint32_t pos=0;pos<n000;pos++)
                    mean_abs_coeff[pos]+=fabs(temp_element1[pos]);
                // 量化
                if (if_zero_freq_expanding)
                    quantize_and_overwrite_array2(temp_element1.data(),element_f_qb,n000,coeff_eb,neg_2_coeff_eb,zero_coeff_bound);
                else
                    quantize_and_overwrite_array1(temp_element1.data(),element_f_qb,n000,coeff_eb,neg_2_coeff_eb);
                // 反变换
                inverse_transform_regular(temp_element1,temp_element2,n0,n00,n000,neg_L0);
                quantize_err_array(temp_element2.data(),element,element_err_qb,n000,abs_eb,neg_2_eb);
            }
            for (uint32_t pos=0;pos<n000;pos++)
                mean_abs_coeff[pos]/=(end_index-begin_index);
            FILE *mean_abs_coeff_file=fopen((data_path+"."+field+".mean_abs_coeff.block"+to_string(block_id)+"_DLT").c_str(),"wb");
            fwrite(mean_abs_coeff.data(),sizeof(T),n000,mean_abs_coeff_file);
            fclose(mean_abs_coeff_file);
            // 无损压缩
            vector<unsigned char> huffman_bs_vec(2*n000*(end_index-begin_index)*sizeof(int64_t),0);
            // {
            //     unsigned char *huffman_bs=huffman_bs_vec.data();
            //     if (if_zero_map)
            //         lossless_compress_coeff_qb_with_zero_map(begin_index,end_index,coeff_qb,huffman_bs,n000);
            //     else
            //         lossless_compress_coeff_qb(begin_index,end_index,coeff_qb,huffman_bs,n000);
            //     uint64_t huffman_bs_len=huffman_bs-huffman_bs_vec.data();
            //     uint64_t max_zstd_bs_len=ZSTD_compressBound(huffman_bs_len);
            //     block.zstd_bs.resize(max_zstd_bs_len);
            //     uint64_t zstd_bs_len=ZSTD_compress(block.zstd_bs.data(),max_zstd_bs_len,huffman_bs_vec.data(),huffman_bs_len,3);
            //     coeff_br[block_id]=zstd_bs_len;
            //     fill(huffman_bs_vec.begin(),huffman_bs_vec.end(),0);
            //     block.zstd_bs.clear();
            // }
            // {
            //     unsigned char *huffman_bs=huffman_bs_vec.data();
            //     huffman_encoder_class<int64_t> err_encoder;
            //     err_encoder.preprocess_encode(err_qb+begin_index*n000,n000*(end_index-begin_index));
            //     err_encoder.save(huffman_bs);
            //     err_encoder.encode(err_qb+begin_index*n000,n000*(end_index-begin_index),huffman_bs);
            //     uint64_t huffman_bs_len=huffman_bs-huffman_bs_vec.data();
            //     uint64_t max_zstd_bs_len=ZSTD_compressBound(huffman_bs_len);
            //     block.zstd_bs.resize(max_zstd_bs_len);
            //     uint64_t zstd_bs_len=ZSTD_compress(block.zstd_bs.data(),max_zstd_bs_len,huffman_bs_vec.data(),huffman_bs_len,3);
            //     err_br[block_id]=zstd_bs_len;
            //     fill(huffman_bs_vec.begin(),huffman_bs_vec.end(),0);
            //     block.zstd_bs.clear();
            // }
            unsigned char *huffman_bs=huffman_bs_vec.data();
            if (if_zero_map)
                lossless_compress_coeff_qb_with_zero_map(begin_index,end_index,coeff_qb,huffman_bs,n000);
            else
                lossless_compress_coeff_qb(begin_index,end_index,coeff_qb,huffman_bs,n000);
            huffman_encoder_class<int64_t> err_encoder;
            err_encoder.preprocess_encode(err_qb+begin_index*n000,n000*(end_index-begin_index));
            err_encoder.save(huffman_bs);
            err_encoder.encode(err_qb+begin_index*n000,n000*(end_index-begin_index),huffman_bs);
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
        printf("\033[91mField= %s, compression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),element_num_to_process*n000*header.word_size/((end_time-start_time)*1024*1024*max_cpus));
        // printf("coeff_br= %f\n",accumulate(coeff_br.begin(),coeff_br.end(),(double)0)*8/element_num_to_process/n000);
        // printf("err_br= %f\n",accumulate(err_br.begin(),err_br.end(),(double)0)*8/element_num_to_process/n000);
    }

    template<typename T>
    inline void param_class<T>::forward_transform_regular(vector<T> &temp_element1,vector<T> &temp_element2,T *element,uint32_t n0,uint32_t n00,uint32_t n000,vector<T> &L0)
    {
        memset(temp_element1.data(),0,sizeof(T)*n000);
        for (uint32_t an0=0;an0<n00;an0+=n0)
            for (uint32_t i=0;i<n0;i++)
                for (uint32_t bn0pc=0;bn0pc<n00;bn0pc++)
                    temp_element1[i*n00+bn0pc]+=element[an0*n0+bn0pc]*L0[an0+i];
        memset(temp_element2.data(),0,sizeof(T)*n000);
        for (uint32_t in00=0;in00<n000;in00+=n00)
            for (uint32_t bn0=0;bn0<n00;bn0+=n0)
                for (uint32_t j=0;j<n0;j++)
                    for (uint32_t c=0;c<n0;c++)
                        temp_element2[in00+j*n0+c]+=temp_element1[in00+bn0+c]*L0[bn0+j];
        memset(temp_element1.data(),0,sizeof(T)*n000);
        for (uint32_t in00pjn0=0;in00pjn0<n000;in00pjn0+=n0)
            for (uint32_t c=0;c<n0;c++)
                for (uint32_t k=0;k<n0;k++)
                    temp_element1[in00pjn0+k]+=temp_element2[in00pjn0+c]*L0[c*n0+k];
    }

    template<typename T>
    inline void param_class<T>::inverse_transform_regular(vector<T> &temp_element1,vector<T> &temp_element2,uint32_t n0,uint32_t n00,uint32_t n000,vector<T> &neg_L0)
    {
        memset(temp_element2.data(),0,sizeof(T)*n000);
        for (uint32_t in00pjn0=0;in00pjn0<n000;in00pjn0+=n0)
            for (uint32_t k=0;k<n0;k++)
                for (uint32_t c=0;c<n0;c++)
                    temp_element2[in00pjn0+c]+=temp_element1[in00pjn0+k]*neg_L0[k*n0+c];
        memset(temp_element1.data(),0,sizeof(T)*n000);
        for (uint32_t in00=0;in00<n000;in00+=n00)
            for (uint32_t jn0=0;jn0<n00;jn0+=n0)
                for (uint32_t b=0;b<n0;b++)
                    for (uint32_t c=0;c<n0;c++)
                        temp_element1[in00+b*n0+c]+=temp_element2[in00+jn0+c]*neg_L0[jn0+b];
        memset(temp_element2.data(),0,sizeof(T)*n000);
        for (uint32_t in0=0;in0<n00;in0+=n0)
            for (uint32_t a=0;a<n0;a++)
                for (uint32_t bn0pc=0;bn0pc<n00;bn0pc++)
                    temp_element2[a*n00+bn0pc]+=temp_element1[in0*n0+bn0pc]*neg_L0[in0+a];
    }

    template<typename T>
    inline void param_class<T>::lossless_compress_coeff_qb_with_zero_map(uint64_t begin_index,uint64_t end_index,int64_t *coeff_qb,unsigned char *&huffman_bs,uint32_t n000)
    {
        bool coeff_qb_all_zero=true;//coeff_qb全0
        bool coeff_qb_none_zero=true;//所有沿着element_id维度的柱子上至少有一个非0
        uint32_t all_zero_num=0;//沿着element_id的维度上全0的柱子数量
        vector<bool> coeff_qb_zero_map(n000,true);//沿着element_id的维度上全0的柱子的标记
        vector<uint32_t> coeff_qb_loose_pos(n000,0);//从紧凑存储时非全0的柱子位置到原本的位置的映射
        for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
            for (uint32_t pos=0;pos<n000;pos++)
                if (coeff_qb[element_id*n000+pos]!=0)
                    coeff_qb_zero_map[pos]=false;
        for (uint32_t pos=0;pos<n000;pos++)
        {
            coeff_qb_loose_pos[pos-all_zero_num]=pos;
            if (coeff_qb_zero_map[pos]==true)
            {
                coeff_qb_none_zero=false;
                all_zero_num++;
            }
            else
                coeff_qb_all_zero=false;
        }
        
        if (coeff_qb_all_zero) // 所有柱子全0
        {
            huffman_bs[0]=0;
            huffman_bs+=1;
        }
        else if (coeff_qb_none_zero) // 没有柱子全0
        {
            huffman_bs[0]=0x80;
            huffman_bs+=1;
            if (if_coefficient_transpose)
            {
                int64_t *temp_coeff_qb=new int64_t[n000*(end_index-begin_index)];
                for (uint32_t pos=0;pos<n000;pos++)
                    for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
                        temp_coeff_qb[pos*(end_index-begin_index)+(element_id-begin_index)]=coeff_qb[element_id*n000+pos];
                huffman_encoder_class<int64_t> f_encoder;
                f_encoder.preprocess_encode(temp_coeff_qb,n000*(end_index-begin_index));
                f_encoder.save(huffman_bs);
                f_encoder.encode(temp_coeff_qb,n000*(end_index-begin_index),huffman_bs);
                delete[] temp_coeff_qb;
            }
            else
            {
                huffman_encoder_class<int64_t> f_encoder;
                f_encoder.preprocess_encode(coeff_qb+begin_index*n000,n000*(end_index-begin_index));
                f_encoder.save(huffman_bs);
                f_encoder.encode(coeff_qb+begin_index*n000,n000*(end_index-begin_index),huffman_bs);
            }
        }
        else // 存在全0柱子也存在非全0柱子
        {
            huffman_bs[0]=0x40;
            uint32_t used_bits=2;
            for (uint32_t pos=0;pos<n000;pos++)
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
            
            if (if_coefficient_transpose)
            {
                int64_t *temp_coeff_qb=new int64_t[(n000-all_zero_num)*(end_index-begin_index)];
                for (uint32_t pos=0;pos<n000-all_zero_num;pos++)
                {
                    uint32_t loose_pos=coeff_qb_loose_pos[pos];
                    for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
                        temp_coeff_qb[pos*(end_index-begin_index)+(element_id-begin_index)]=coeff_qb[element_id*n000+loose_pos];
                }
                huffman_encoder_class<int64_t> f_encoder;
                f_encoder.preprocess_encode(temp_coeff_qb,(n000-all_zero_num)*(end_index-begin_index));
                f_encoder.save(huffman_bs);
                f_encoder.encode(temp_coeff_qb,(n000-all_zero_num)*(end_index-begin_index),huffman_bs);
                delete[] temp_coeff_qb;
            }
            else
            {
                for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
                {
                    auto original_element_f_qb=coeff_qb+element_id*n000;
                    auto compact_element_f_qb=coeff_qb+begin_index*n000+(element_id-begin_index)*(n000-all_zero_num);
                    for (uint32_t pos=0;pos<(n000-all_zero_num);pos++)
                        compact_element_f_qb[pos]=original_element_f_qb[coeff_qb_loose_pos[pos]];
                }
                huffman_encoder_class<int64_t> f_encoder;
                f_encoder.preprocess_encode(coeff_qb+begin_index*n000,(n000-all_zero_num)*(end_index-begin_index));
                f_encoder.save(huffman_bs);
                f_encoder.encode(coeff_qb+begin_index*n000,(n000-all_zero_num)*(end_index-begin_index),huffman_bs);
            }
        }   
    }

    template<typename T>
    inline void param_class<T>::lossless_compress_coeff_qb(uint64_t begin_index,uint64_t end_index,int64_t *coeff_qb,unsigned char *&huffman_bs,uint32_t n000)
    {
        if (if_coefficient_transpose)
        {
            int64_t *temp_coeff_qb=new int64_t[n000*(end_index-begin_index)];
            for (uint32_t pos=0;pos<n000;pos++)
                for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
                    temp_coeff_qb[pos*(end_index-begin_index)+(element_id-begin_index)]=coeff_qb[element_id*n000+pos];
            huffman_encoder_class<int64_t> f_encoder;
            f_encoder.preprocess_encode(temp_coeff_qb,n000*(end_index-begin_index));
            f_encoder.save(huffman_bs);
            f_encoder.encode(temp_coeff_qb,n000*(end_index-begin_index),huffman_bs);
            delete[] temp_coeff_qb;
        }
        else
        {
            huffman_encoder_class<int64_t> f_encoder;
            f_encoder.preprocess_encode(coeff_qb+begin_index*n000,n000*(end_index-begin_index));
            f_encoder.save(huffman_bs);
            f_encoder.encode(coeff_qb+begin_index*n000,n000*(end_index-begin_index),huffman_bs);
        }
    }
}
