#pragma once

# XXX 本文件归档已移除的partition级zero map解压实现，不参与编译。
# XXX 对应zero_map_compress.tpp；实验收益过小，详情见该文件开头的说明。

namespace mrcz
{
    template<typename T>
    inline void param_class<T>::lossless_decompress_coeff_qb_with_zero_map(uint64_t begin_index,uint64_t end_index,int64_t *coeff_qb,const unsigned char *&huffman_bs,uint32_t n000)
    {
        uint32_t all_zero_num=0;
        vector<bool> coeff_qb_zero_map(n000,false);
        vector<uint32_t> coeff_qb_loose_pos(n000,0),coeff_qb_compact_pos(n000,0);
        if ((huffman_bs[0]&0xC0)==0)
        {
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
                for (uint32_t pos=0;pos<n000;pos++)
                {
                    uint32_t compact_pos=coeff_qb_compact_pos[pos];
                    if (coeff_qb_zero_map[pos]==false)
                        for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
                            coeff_qb[element_id*n000+pos]=temp_coeff_qb[compact_pos*(end_index-begin_index)+(element_id-begin_index)];
                    else
                        for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
                            coeff_qb[element_id*n000+pos]=0;
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
