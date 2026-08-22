#pragma once

#include "baseline_compress.hpp"
#include "mrcz_core/mrcz.hpp"

namespace mrcz
{
    template<typename T>
    void compress(param_class<T> &params)
    {
        if (params.method.empty() or params.method=="MRCZ-3D-GLL" or params.method=="MRCZ-3D-GLL-experimental" or
            params.method=="SZ3-3D-GLL" or params.method=="SZ3-3D-GLL-experimental" or params.method=="QoZ-3D-GLL")
            core::compress(params);
        else
            baseline::compress(params);
    }
}
