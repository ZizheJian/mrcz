#pragma once

namespace mrcz
{
    inline bool is_mrcz_coefficient_delta_method(string method)
    {
        return method=="MRCZ-3D-GLL-coeff-delta-constant" or method=="MRCZ-3D-GLL-coeff-delta-linear" or method=="MRCZ-3D-GLL-coeff-delta-all";
    }

    template<typename T>
    inline int param_class<T>::coefficient_delta_max_order()
    {
        if (method=="MRCZ-3D-GLL-coeff-delta-constant")
            return 0;
        if (method=="MRCZ-3D-GLL-coeff-delta-linear")
            return 1;
        if (method=="MRCZ-3D-GLL-coeff-delta-all")
            return numeric_limits<int>::max();
        return -1;
    }

    template<typename T>
    inline void param_class<T>::apply_coefficient_delta(int64_t *coefficient_major_qb,uint64_t block_element_num)
    {
        int max_order=coefficient_delta_max_order();
        if (max_order<0 or block_element_num<2)
            return;
        uint32_t n0=header.c_element_shape[0],n1=header.c_element_shape[1],n2=header.c_element_shape[2];
        for (uint32_t i=0;i<n0;i++)
            for (uint32_t j=0;j<n1;j++)
                for (uint32_t k=0;k<n2;k++)
                {
                    if (i+j+k>static_cast<uint32_t>(max_order))
                        continue;
                    uint32_t pos=i*n1*n2+j*n2+k;
                    int64_t *mode_qb=coefficient_major_qb+pos*block_element_num;
                    for (uint64_t element_id=block_element_num-1;element_id>0;element_id--)
                    {
                        uint64_t current=bit_cast<uint64_t>(mode_qb[element_id]);
                        uint64_t previous=bit_cast<uint64_t>(mode_qb[element_id-1]);
                        mode_qb[element_id]=bit_cast<int64_t>(current-previous);
                    }
                }
    }

    template<typename T>
    inline void param_class<T>::undo_coefficient_delta(int64_t *coefficient_major_qb,uint64_t block_element_num)
    {
        int max_order=coefficient_delta_max_order();
        if (max_order<0 or block_element_num<2)
            return;
        uint32_t n0=header.c_element_shape[0],n1=header.c_element_shape[1],n2=header.c_element_shape[2];
        for (uint32_t i=0;i<n0;i++)
            for (uint32_t j=0;j<n1;j++)
                for (uint32_t k=0;k<n2;k++)
                {
                    if (i+j+k>static_cast<uint32_t>(max_order))
                        continue;
                    uint32_t pos=i*n1*n2+j*n2+k;
                    int64_t *mode_qb=coefficient_major_qb+pos*block_element_num;
                    for (uint64_t element_id=1;element_id<block_element_num;element_id++)
                    {
                        uint64_t delta=bit_cast<uint64_t>(mode_qb[element_id]);
                        uint64_t previous=bit_cast<uint64_t>(mode_qb[element_id-1]);
                        mode_qb[element_id]=bit_cast<int64_t>(delta+previous);
                    }
                }
    }
}
