#pragma once

#include "baseline_compress.hpp"
#include "mrcz_core/mrcz.hpp"

namespace mrcz
{
    inline void reject_old_compress_method_name(string method)
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
    void compress(param_class<T> &params)
    {
        reject_old_compress_method_name(params.method);
        if (params.method.empty() or params.method=="MRCZ-3D-GLL" or params.method=="MRCZ-3D-GLL-experimental" or
            params.method=="SZ3-3D-GLL" or params.method=="QoZ-3D-GLL")
            core::compress(params);
        else
            baseline::compress(params);
    }
}
