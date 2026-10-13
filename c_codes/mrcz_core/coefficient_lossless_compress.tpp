#pragma once

namespace mrcz
{
    template<typename T>
    template<bool ablation>
    inline void param_class<T>::lossless_compress_coeff_qb(uint64_t begin_index,uint64_t end_index,int32_t *coeff_qb,unsigned char *&huffman_bs,uint32_t n000)
    {
        bool coefficient_transpose=true,tiled_coefficient_transpose=true;
        if constexpr(ablation)
        {
            coefficient_transpose=if_coefficient_transpose;
            tiled_coefficient_transpose=if_tiled_coefficient_transpose;
        }
        if (coefficient_transpose)
        {
            uint64_t block_element_num=end_index-begin_index;
            int32_t *temp_coeff_qb=new int32_t[n000*block_element_num];
            if (tiled_coefficient_transpose)
            {
                constexpr uint32_t tile_size=32;
                for (uint32_t pos_begin=0;pos_begin<n000;pos_begin+=tile_size)
                    for (uint64_t element_begin=begin_index;element_begin<end_index;element_begin+=tile_size)
                        for (uint32_t pos=pos_begin;pos<min<uint32_t>(pos_begin+tile_size,n000);pos++)
                            for (uint64_t element_id=element_begin;element_id<min<uint64_t>(element_begin+tile_size,end_index);element_id++)
                                temp_coeff_qb[pos*block_element_num+(element_id-begin_index)]=coeff_qb[element_id*n000+pos];
            }
            else
                for (uint32_t pos=0;pos<n000;pos++)
                    for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
                        temp_coeff_qb[pos*block_element_num+(element_id-begin_index)]=coeff_qb[element_id*n000+pos];
            fast_huffman_encoder_class<int32_t> coeff_encoder;
            coeff_encoder.preprocess_encode(temp_coeff_qb,n000*block_element_num);
            coeff_encoder.save(huffman_bs);
            coeff_encoder.encode(temp_coeff_qb,n000*block_element_num,huffman_bs);
            delete[] temp_coeff_qb;
        }
        else
        {
            fast_huffman_encoder_class<int32_t> coeff_encoder;
            coeff_encoder.preprocess_encode(coeff_qb+begin_index*n000,n000*(end_index-begin_index));
            coeff_encoder.save(huffman_bs);
            coeff_encoder.encode(coeff_qb+begin_index*n000,n000*(end_index-begin_index),huffman_bs);
        }
    }
}
