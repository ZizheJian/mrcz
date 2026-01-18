#pragma once

namespace jzzz
{
    template<typename T>
    void param_class<T>::quantize_and_overwrite(T &data,int64_t &qb,T eb,T neg_2_eb)
    {
        qb=static_cast<int64_t>(round(data*neg_2_eb));
        data=qb*2*eb;
    }

    template<typename T>
    void param_class<T>::quantize_and_overwrite_err(T &pred,T &truth,int64_t &qb,T eb,T neg_2_eb)
    {
        qb=static_cast<int64_t>(round((pred-truth)*neg_2_eb));
        truth=pred-qb*2*eb;
    }

    template<typename T>
    void param_class<T>::quantize_err_array(T *pred,T *truth,int64_t *qb,uint32_t size,T eb,T neg_2_eb)
    {
        for (uint32_t i=0;i<size;i++)
            qb[i]=static_cast<int64_t>(round((pred[i]-truth[i])*neg_2_eb));
    }

    template<typename T>
    void param_class<T>::quantize_and_overwrite_array2(T *data,int64_t *qb,uint32_t size,T eb,T neg_2_eb,T zero_freq_bound)
    {
        for (uint32_t i=0;i<size;i++)
        {
            if (-zero_freq_bound<=data[i] && data[i]<=zero_freq_bound)
                qb[i]=0;
            else
                qb[i]=static_cast<int64_t>(round(data[i]*neg_2_eb));
            data[i]=qb[i]*2*eb;
        }
    }
}