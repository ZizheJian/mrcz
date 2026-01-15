#ifndef __QUANTIZE_CPP__
#define __QUANTIZE_CPP__

#include <cmath>
#include "jzzz_param.hpp"

namespace jzzz
{
    void param_class::quantize_and_overwrite_number(float &data,int64_t &qb,float eb,float neg_2_eb)
    {
        qb=static_cast<int64_t>(round(data*neg_2_eb));
        data=qb*2*eb;
    }

    void param_class::quantize_err(float *pred,float *truth,int64_t *qb,uint32_t size,float eb,float neg_2_eb)
    {
        for (uint32_t i=0;i<size;i++)
            qb[i]=static_cast<int64_t>(round((pred[i]-truth[i])*neg_2_eb));
    }

    void param_class::quantize_and_overwrite2(float *data,int64_t *qb,uint32_t size,float eb,float neg_2_eb,float zero_freq_bound)
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

#endif