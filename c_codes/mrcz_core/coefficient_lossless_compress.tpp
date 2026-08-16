#pragma once

namespace mrcz
{
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

