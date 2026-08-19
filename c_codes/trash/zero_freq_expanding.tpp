#pragma once

// XXX 废弃的zero-frequency expansion实验代码。
// 该方法把量化为0的系数区间从普通的[-eb,eb]扩大到
// [-zero_coeff_qb_bound*eb,zero_coeff_qb_bound*eb]，原默认倍数为3。
// 实验中，它在相同error bound下可以降低bitrate，但也会明显降低PSNR。
// 在mixlay0.f00001和jicf0.f00000上按相同PSNR比较时，压缩性能均比基础MRCZ更差，
// 因此从mrcz-experimental中移除。这里仅保留代码，方便以后追溯，不参与构建。

// 原压缩端调用逻辑：
// T zero_coeff_bound=coeff_eb*zero_coeff_qb_bound;
// if (method=="mrcz-experimental" and if_zero_freq_expanding)
//     quantize_and_overwrite_array2(data,qb,size,coeff_eb,neg_2_coeff_eb,zero_coeff_bound);
// else
//     quantize_and_overwrite_array1(data,qb,size,coeff_eb,neg_2_coeff_eb);

namespace mrcz
{
    template<typename T>
    inline void param_class<T>::quantize_and_overwrite_array2(T *data,int64_t *qb,uint32_t size,T eb,T neg_2_eb,T zero_coeff_bound)
    {
        for (uint32_t i=0;i<size;i++)
        {
            if (-zero_coeff_bound<=data[i] and data[i]<=zero_coeff_bound)
                qb[i]=0;
            else
                qb[i]=static_cast<int64_t>(round(data[i]*neg_2_eb));
            data[i]=qb[i]*2*eb;
        }
    }
}
