#pragma once

namespace mrcz
{
    template<typename T>
    inline void param_class<T>::quantize_and_overwrite(T &data,int64_t &qb,T eb,T neg_2_eb)
    {
        qb=static_cast<int64_t>(round(data*neg_2_eb));
        data=qb*2*eb;
    }

    template<typename T>
    inline void param_class<T>::quantize_and_overwrite_err(T &pred,T &truth,int64_t &qb,T eb,T neg_2_eb)
    {
        qb=static_cast<int64_t>(round((pred-truth)*neg_2_eb));
        truth=pred-qb*2*eb;
    }

    template<typename T>
    inline void param_class<T>::quantize_err_array(T *pred,T *truth,int64_t *qb,uint32_t size,T eb,T neg_2_eb)
    {
        for (uint32_t i=0;i<size;i++)
            qb[i]=static_cast<int64_t>(round((pred[i]-truth[i])*neg_2_eb));
    }

    template<typename T>
    inline void param_class<T>::quantize_and_overwrite_array1(T *data,int64_t *qb,uint32_t size,T eb,T neg_2_eb)
    {
        for (uint32_t i=0;i<size;i++)
        {
            qb[i]=static_cast<int64_t>(round(data[i]*neg_2_eb));
            data[i]=qb[i]*2*eb;
        }
    }

}
