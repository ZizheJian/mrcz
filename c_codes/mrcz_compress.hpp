#pragma once

#include "baseline_compress.hpp"
#include "mrcz_core/compress.hpp"

namespace mrcz
{
    template<typename T>
    void compress(param_class<T> &params)
    {
        if (params.method.empty() or params.method=="mrcz" or params.method=="sz3-hohs")
            core::compress(params);
        else
            baseline::compress(params);
    }
}
