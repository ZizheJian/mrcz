#pragma once

#include "baseline_decompress.hpp"
#include "mrcz_core/decompress.hpp"

namespace mrcz
{
    template<typename T>
    void decompress(param_class<T> &params)
    {
        if (params.method.empty() or params.method=="mrcz" or params.method=="sz3-hohs")
            core::decompress(params);
        else
            baseline::decompress(params);
    }
}
