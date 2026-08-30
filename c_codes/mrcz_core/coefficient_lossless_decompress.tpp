#pragma once

namespace mrcz
{
    template<typename T>
    inline void param_class<T>::lossless_decompress_coeff_qb(uint64_t begin_index,uint64_t end_index,int64_t *coeff_qb,const unsigned char *&huffman_bs,uint32_t n000)
    {
        if (if_coefficient_transpose)
        {
            uint64_t block_element_num=end_index-begin_index;
            int64_t *temp_coeff_qb=new int64_t[n000*block_element_num];
            huffman_encoder_class<int64_t> coeff_encoder;
            coeff_encoder.load(huffman_bs);
            coeff_encoder.decode(huffman_bs,temp_coeff_qb,n000*block_element_num);
            if (if_tiled_coefficient_transpose)
            {
                constexpr uint32_t tile_size=32;
                for (uint32_t pos_begin=0;pos_begin<n000;pos_begin+=tile_size)
                    for (uint64_t element_begin=begin_index;element_begin<end_index;element_begin+=tile_size)
                        for (uint32_t pos=pos_begin;pos<min<uint32_t>(pos_begin+tile_size,n000);pos++)
                            for (uint64_t element_id=element_begin;element_id<min<uint64_t>(element_begin+tile_size,end_index);element_id++)
                                coeff_qb[element_id*n000+pos]=temp_coeff_qb[pos*block_element_num+(element_id-begin_index)];
            }
            else
                for (uint32_t pos=0;pos<n000;pos++)
                    for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
                        coeff_qb[element_id*n000+pos]=temp_coeff_qb[pos*block_element_num+(element_id-begin_index)];
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
