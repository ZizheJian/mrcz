#pragma once

#include "baseline_compress.hpp"
#include "mrcz_core/mrcz.hpp"

namespace mrcz
{
    inline void reject_old_compress_method_name(string method)
    {
        if (method=="mrcz-experimental" or method=="mrcz-3d-gll-experimental")
            throw runtime_error("method "+method+"已重命名为MRCZ-prediction-diagnostics；该方法只做预测诊断，不生成bitstream");
    }

    template<typename T>
    void compress(param_class<T> &params)
    {
        reject_old_compress_method_name(params.method);
        if (params.method.empty() or params.method=="MRCZ" or params.method=="MRCZ-element-batched" or params.method=="MRCZ-element-batched-original-transpose" or params.method=="MRCZ-prediction-diagnostics" or
            params.method=="FDWZ" or params.method=="FDWZ-avx2-even" or params.method=="FDWZ-avx2-away" or params.method=="FDWZ-five-trees" or params.method=="FDWZ-predictability" or params.method=="FDWZ-fixed-cubic" or params.method=="FDWZ-surface-first" or params.method=="FDWZ-normal-boundary" or params.method=="FDWZ-index-order" or params.method=="FDWZ-index-midpoint" or params.method=="QoZ-3D-GLL" or
            params.method.rfind("FDWZ-order-",0)==0 or params.method.rfind("FDWZ-reference-EB-",0)==0 or params.method.rfind("FDWZ-sampled-",0)==0 or params.method.rfind("FDWZ-position-bias-",0)==0)
            core::compress(params);
        else
            baseline::compress(params);
    }
}
