#pragma once

#include "baseline_compress.hpp"
#include "mrcz_core/mrcz.hpp"

namespace mrcz
{
    template<typename T>
    void compress(param_class<T> &params)
    {
        if (params.method.empty() or params.method=="MRCZ-3D-GLL" or params.method=="MRCZ-3D-GLL-experimental" or
            params.method=="SZ3-3D-GLL" or params.method=="SZ3-3D-GLL-optimized" or params.method=="SZ3-3D-GLL-optimized-avx2-even" or params.method=="SZ3-3D-GLL-optimized-avx2-away" or params.method=="SZ3-3D-GLL-experimental" or params.method=="SZ3-3D-GLL-predictability" or params.method=="SZ3-3D-GLL-cubic" or params.method=="SZ3-3D-GLL-surface-first" or params.method=="SZ3-3D-GLL-normal-boundary" or params.method=="SZ3-3D-GLL-index-order" or params.method=="SZ3-3D-GLL-index-midpoint" or params.method=="QoZ-3D-GLL" or
            params.method.rfind("SZ3-3D-GLL-order-",0)==0 or params.method.rfind("SZ3-3D-GLL-reference-EB-",0)==0 or params.method.rfind("SZ3-3D-GLL-sampled-",0)==0 or params.method.rfind("SZ3-3D-GLL-position-bias-",0)==0)
            core::compress(params);
        else
            baseline::compress(params);
    }
}
