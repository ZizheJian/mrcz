#pragma once

#include "baseline_compress.hpp"
#include "mrcz_core/mrcz.hpp"

namespace mrcz
{
    inline void reject_old_compress_method_name(string method)
    {
        if (method=="mrcz") throw runtime_error("method mrcz已重命名为mrcz-3d-gll");
        if (method=="mrcz-experimental") throw runtime_error("method mrcz-experimental已重命名为mrcz-3d-gll-experimental");
        if (method=="sz3-hohs") throw runtime_error("method sz3-hohs已重命名为sz3-3d-gll");
        if (method=="sz3-1") throw runtime_error("method sz3-1已重命名为sz3-3d-uniform");
        if (method=="sz2-0") throw runtime_error("method sz2-0已重命名为sz2-3d-uniform");
        if (method=="sz2-1") throw runtime_error("method sz2-1已重命名为sz2-3d-gll");
        if (method=="zfp") throw runtime_error("method zfp已重命名为zfp-4d-uniform");
        if (method=="sperr") throw runtime_error("method sperr已重命名为sperr-3d-uniform");
    }

    template<typename T>
    void compress(param_class<T> &params)
    {
        reject_old_compress_method_name(params.method);
        if (params.method.empty() or params.method=="mrcz-3d-gll" or params.method=="mrcz-3d-gll-experimental" or params.method=="sz3-3d-gll")
            core::compress(params);
        else
            baseline::compress(params);
    }
}
