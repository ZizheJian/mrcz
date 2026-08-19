#pragma once

namespace mrcz
{
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
