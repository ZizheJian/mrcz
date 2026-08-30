#pragma once

#include "baseline_decompress.hpp"
#include "mrcz_compress.hpp"

namespace mrcz
{
    template<typename T>
    void decompress(param_class<T> &params)
    {
        if (is_core_method(params.method))
            core::decompress(params);
        else if (is_baseline_method(params.method))
            baseline::decompress(params);
        else
            throw runtime_error("未知的解压方法: "+params.method);
    }
}
