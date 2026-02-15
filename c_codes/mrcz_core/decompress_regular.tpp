#pragma once

namespace mrcz
{
    template<typename T>
    void param_class<T>::decompress_regular(string field)
    {
        uint32_t n0=header.element_shape[0];
        uint32_t n00=n0*n0;
        uint32_t n000=n00*n0;
        vector<T> L0=L[0];
        T abs_eb=this->abs_eb[field];
        T coeff_eb=abs_eb*coeff_eb_scale;
        int64_t *coeff_qb=new int64_t[element_num_to_process*n000];
        int64_t *err_qb=new int64_t[element_num_to_process*n000];
        decompressed_data.resize(element_num_to_process*n000);
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
                lossless_decompress_coeff_qb_with_zero_map(begin_index,end_index,coeff_qb,huffman_bs,n000);
            else
                lossless_decompress_coeff_qb(begin_index,end_index,coeff_qb,huffman_bs,n000);
            huffman_encoder_class<int64_t> err_encoder;
            err_encoder.load(huffman_bs);
            err_encoder.decode(huffman_bs,err_qb+n000*begin_index,n000*(end_index-begin_index));
            
            // 然后重构剩余的数据
            vector<T> temp_element1(n000,0),temp_element2(n000,0);
            for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
            {
                T *element=decompressed_data.data()+element_id*n000;
                int64_t *element_coeff_qb=coeff_qb+element_id*n000;
                int64_t *element_err_qb=err_qb+element_id*n000;
                // 重构coeff
                for (uint32_t i=0;i<n000;i++)
                    temp_element1[i]=element_coeff_qb[i]*2*coeff_eb;
                // 反变换
                inverse_transform_regular(temp_element1,temp_element2,n0,n00,n000,L0);
                // 补偿误差
                for (uint32_t j=0;j<n000;j++)
                    element[j]+=(temp_element2[j]-element_err_qb[j]*2*abs_eb);
            }
        }
        double end_time=omp_get_wtime();
        delete[] err_qb;
        delete[] coeff_qb;
        printf("Decompression time= %f seconds\n",end_time-start_time);
        printf("\033[91mBandwidth per cpu= %f MB/s\033[0m\n",element_num_to_process*n000*header.word_size/((end_time-start_time)*1024*1024*max_cpus));
    }

    template<typename T>
    inline void param_class<T>::lossless_decompress_coeff_qb_with_zero_map(uint64_t begin_index,uint64_t end_index,int64_t *coeff_qb,const unsigned char *&huffman_bs,uint32_t n000)
    {
        bool coeff_qb_all_zero=false;
        uint32_t all_zero_num=0;
        vector<bool> coeff_qb_zero_map(n000,false);
        vector<uint32_t> coeff_qb_loose_pos(n000,0),coeff_qb_compact_pos(n000,0);
        if ((huffman_bs[0]&0xC0)==0)
        {
            coeff_qb_all_zero=true;
            coeff_qb_zero_map.assign(n000,true);
            huffman_bs+=1;
            memset(coeff_qb+begin_index*n000,0,sizeof(int64_t)*n000*(end_index-begin_index));
        }
        else if ((huffman_bs[0]&0xC0)==0x80)
        {
            huffman_bs+=1;
            if (if_coefficient_transpose)
            {
                int64_t *temp_coeff_qb=new int64_t[n000*(end_index-begin_index)];
                huffman_encoder_class<int64_t> coeff_encoder;
                coeff_encoder.load(huffman_bs);
                coeff_encoder.decode(huffman_bs,temp_coeff_qb,n000*(end_index-begin_index));
                for (uint32_t pos=0;pos<n000;pos++)
                    for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
                        coeff_qb[element_id*n000+pos]=temp_coeff_qb[pos*(end_index-begin_index)+(element_id-begin_index)];
                delete[] temp_coeff_qb;
            }
            else
            {
                huffman_encoder_class<int64_t> coeff_encoder;
                coeff_encoder.load(huffman_bs);
                coeff_encoder.decode(huffman_bs,coeff_qb+n000*begin_index,n000*(end_index-begin_index));
            }
        }
        else
        {
            uint32_t used_bits=2;
            for (uint32_t pos=0;pos<n000;pos++)
            {
                if ((huffman_bs[0]>>(7-used_bits))&0x01)
                    coeff_qb_zero_map[pos]=true;
                used_bits++;
                if (used_bits==8)
                {
                    huffman_bs++;
                    used_bits=0;
                }
            }
            if (used_bits!=0)
                huffman_bs++;
            
            for (uint32_t pos=0;pos<n000;pos++)
            {
                coeff_qb_loose_pos[pos-all_zero_num]=pos;
                coeff_qb_compact_pos[pos]=pos-all_zero_num;
                if (coeff_qb_zero_map[pos])
                    all_zero_num++;
            }
            if (if_coefficient_transpose)
            {
                int64_t *temp_coeff_qb=new int64_t[(n000-all_zero_num)*(end_index-begin_index)];
                huffman_encoder_class<int64_t> coeff_encoder;
                coeff_encoder.load(huffman_bs);
                coeff_encoder.decode(huffman_bs,temp_coeff_qb,(n000-all_zero_num)*(end_index-begin_index));
                for (uint32_t pos=0;pos<n000-all_zero_num;pos++)
                {
                    uint32_t loose_pos=coeff_qb_loose_pos[pos];
                    for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
                        coeff_qb[element_id*n000+loose_pos]=temp_coeff_qb[pos*(end_index-begin_index)+(element_id-begin_index)];
                }
                delete[] temp_coeff_qb;
            }
            else
            {
                huffman_encoder_class<int64_t> coeff_encoder;
                coeff_encoder.load(huffman_bs);
                coeff_encoder.decode(huffman_bs,coeff_qb+n000*begin_index,(n000-all_zero_num)*(end_index-begin_index));
                for (uint64_t neg_element_id=begin_index;neg_element_id<end_index;neg_element_id++)
                {
                    uint64_t element_id=begin_index+end_index-1-neg_element_id;
                    auto compact_element_coeff_qb=coeff_qb+begin_index*n000+(element_id-begin_index)*(n000-all_zero_num);
                    auto original_element_coeff_qb=coeff_qb+element_id*n000;
                    for (uint32_t neg_pos=0;neg_pos<n000;neg_pos++)
                    {
                        uint32_t pos=n000-1-neg_pos;
                        if (coeff_qb_zero_map[pos]==false)
                            original_element_coeff_qb[pos]=compact_element_coeff_qb[coeff_qb_compact_pos[pos]];
                        else
                            original_element_coeff_qb[pos]=0;
                    }
                }
            }
        }
            
    }

    template<typename T>
    inline void param_class<T>::lossless_decompress_coeff_qb(uint64_t begin_index,uint64_t end_index,int64_t *coeff_qb,const unsigned char *&huffman_bs,uint32_t n000)
    {
        if (if_coefficient_transpose)
        {
            int64_t *temp_coeff_qb=new int64_t[(n000)*(end_index-begin_index)];
            huffman_encoder_class<int64_t> coeff_encoder;
            coeff_encoder.load(huffman_bs);
            coeff_encoder.decode(huffman_bs,temp_coeff_qb,n000*(end_index-begin_index));
            for (uint32_t pos=0;pos<n000;pos++)
                for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
                    coeff_qb[element_id*n000+pos]=temp_coeff_qb[pos*(end_index-begin_index)+(element_id-begin_index)];
            delete[] temp_coeff_qb;
        }
        else
        {
            huffman_encoder_class<int64_t> coeff_encoder;
            coeff_encoder.load(huffman_bs);
            coeff_encoder.decode(huffman_bs,coeff_qb+n000*begin_index,n000*(end_index-begin_index));
        }
    }
}