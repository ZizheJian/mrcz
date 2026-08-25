#pragma once

#include "baseline_decompress.hpp"
#include "mrcz_core/mrcz.hpp"

namespace mrcz
{
    inline void reject_old_decompress_method_name(string method)
    {
        if (method=="mrcz") throw runtime_error("method mrcz已重命名为MRCZ-3D-GLL");
        if (method=="mrcz-experimental") throw runtime_error("method mrcz-experimental已重命名为MRCZ-3D-GLL-experimental");
        if (method=="sz3-hohs") throw runtime_error("method sz3-hohs已重命名为SZ3-3D-GLL");
        if (method=="sz3-1") throw runtime_error("method sz3-1已重命名为SZ3-3D-uniform");
        if (method=="sz2-0") throw runtime_error("method sz2-0已重命名为SZ2-3D-uniform");
        if (method=="sz2-1") throw runtime_error("method sz2-1已重命名为SZ2-3D-GLL");
        if (method=="zfp") throw runtime_error("method zfp已重命名为ZFP-4D-uniform");
        if (method=="sperr") throw runtime_error("method sperr已重命名为SPERR-3D-uniform");
        if (method=="mrcz-3d-gll") throw runtime_error("method mrcz-3d-gll已重命名为MRCZ-3D-GLL");
        if (method=="mrcz-3d-gll-experimental") throw runtime_error("method mrcz-3d-gll-experimental已重命名为MRCZ-3D-GLL-experimental");
        if (method=="sz3-3d-uniform") throw runtime_error("method sz3-3d-uniform已重命名为SZ3-3D-uniform");
        if (method=="sz3-3d-gll") throw runtime_error("method sz3-3d-gll已重命名为SZ3-3D-GLL");
        if (method=="sz3-3d-gll-experimental") throw runtime_error("method sz3-3d-gll-experimental已重命名为SZ3-3D-GLL-experimental");
        if (method=="sz3-3d-gll-surface-first") throw runtime_error("method sz3-3d-gll-surface-first已重命名为SZ3-3D-GLL-surface-first");
        if (method=="sz3-3d-gll-normal-boundary") throw runtime_error("method sz3-3d-gll-normal-boundary已重命名为SZ3-3D-GLL-normal-boundary");
        if (method=="qoz-3d-gll") throw runtime_error("method qoz-3d-gll已重命名为QoZ-3D-GLL");
        if (method=="sz3-1d-uniform") throw runtime_error("method sz3-1d-uniform已重命名为SZ3-1D-uniform");
        if (method=="sz3-4d-uniform") throw runtime_error("method sz3-4d-uniform已重命名为SZ3-4D-uniform");
        if (method=="qoz-4d-uniform") throw runtime_error("method qoz-4d-uniform已重命名为QoZ-4D-uniform");
        if (method=="sz2-3d-uniform") throw runtime_error("method sz2-3d-uniform已重命名为SZ2-3D-uniform");
        if (method=="sz2-3d-gll") throw runtime_error("method sz2-3d-gll已重命名为SZ2-3D-GLL");
        if (method=="zfp-4d-uniform") throw runtime_error("method zfp-4d-uniform已重命名为ZFP-4D-uniform");
        if (method=="sperr-3d-uniform") throw runtime_error("method sperr-3d-uniform已重命名为SPERR-3D-uniform");
    }

    template<typename T>
    void decompress(param_class<T> &params)
    {
        reject_old_decompress_method_name(params.method);
        if (params.method.empty() or params.method=="MRCZ-3D-GLL" or params.method=="MRCZ-3D-GLL-experimental" or
            params.method=="SZ3-3D-GLL" or params.method=="SZ3-3D-GLL-optimized" or params.method=="SZ3-3D-GLL-experimental" or params.method=="SZ3-3D-GLL-cubic" or params.method=="SZ3-3D-GLL-surface-first" or params.method=="SZ3-3D-GLL-normal-boundary" or params.method=="SZ3-3D-GLL-index-order" or params.method=="SZ3-3D-GLL-index-midpoint" or params.method=="QoZ-3D-GLL" or
            params.method.rfind("SZ3-3D-GLL-order-",0)==0 or params.method.rfind("SZ3-3D-GLL-reference-EB-",0)==0 or params.method.rfind("SZ3-3D-GLL-sampled-",0)==0 or params.method.rfind("SZ3-3D-GLL-position-bias-",0)==0)
            core::decompress(params);
        else
            baseline::decompress(params);
    }
}
