#pragma once

#include "baseline_compress.hpp"
#include "mrcz_core/mrcz.hpp"

namespace mrcz
{
    inline void reject_old_compress_method_name(string method)
    {
        if (method=="MRCZ-3D-GLL-experimental" or method=="mrcz-experimental" or method=="mrcz-3d-gll-experimental")
            throw runtime_error("method "+method+"已重命名为MRCZ-prediction-diagnostics；该方法只做预测诊断，不生成bitstream");
    }

    template<typename T>
    void compress(param_class<T> &params)
    {
        reject_old_compress_method_name(params.method);
        if (params.method.empty() or params.method=="MRCZ-3D-GLL" or params.method=="MRCZ-3D-GLL-element-batched" or params.method=="MRCZ-3D-GLL-element-batched-original-transpose" or params.method=="MRCZ-prediction-diagnostics" or
            params.method=="SZ3-3D-GLL" or params.method=="SZ3-3D-GLL-optimized" or params.method=="SZ3-3D-GLL-optimized-avx2-even" or params.method=="SZ3-3D-GLL-optimized-avx2-away" or params.method=="SZ3-3D-GLL-experimental" or params.method=="SZ3-3D-GLL-predictability" or params.method=="SZ3-3D-GLL-cubic" or params.method=="SZ3-3D-GLL-surface-first" or params.method=="SZ3-3D-GLL-normal-boundary" or params.method=="SZ3-3D-GLL-index-order" or params.method=="SZ3-3D-GLL-index-midpoint" or params.method=="QoZ-3D-GLL" or
            params.method.rfind("SZ3-3D-GLL-order-",0)==0 or params.method.rfind("SZ3-3D-GLL-reference-EB-",0)==0 or params.method.rfind("SZ3-3D-GLL-sampled-",0)==0 or params.method.rfind("SZ3-3D-GLL-position-bias-",0)==0)
            core::compress(params);
        else
            baseline::compress(params);
    }
}
