#pragma once

#include "baseline_decompress.hpp"
#include "mrcz_core/mrcz.hpp"

namespace mrcz
{
    inline void reject_old_decompress_method_name(string method)
    {
        if (method=="MRCZ-prediction-diagnostics" or method=="mrcz-experimental" or method=="mrcz-3d-gll-experimental")
            throw runtime_error("MRCZ prediction diagnostics不生成bitstream，不能解压");
        if (method=="mrcz") throw runtime_error("method mrcz已重命名为MRCZ");
        if (method=="sz3-1") throw runtime_error("method sz3-1已重命名为SZ3-3D-uniform");
        if (method=="sz2-0") throw runtime_error("method sz2-0已重命名为SZ2-3D-uniform");
        if (method=="sz2-1") throw runtime_error("method sz2-1已重命名为SZ2-3D-GLL");
        if (method=="zfp") throw runtime_error("method zfp已重命名为ZFP-4D-uniform");
        if (method=="sperr") throw runtime_error("method sperr已重命名为SPERR-3D-uniform");
        if (method=="mrcz-3d-gll") throw runtime_error("method mrcz-3d-gll已重命名为MRCZ");
        if (method=="sz3-3d-uniform") throw runtime_error("method sz3-3d-uniform已重命名为SZ3-3D-uniform");
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
        if (params.method.empty() or params.method=="MRCZ" or params.method=="MRCZ-element-batched" or params.method=="MRCZ-element-batched-original-transpose" or
            params.method=="FDWZ" or params.method=="FDWZ-avx2-even" or params.method=="FDWZ-avx2-away" or params.method=="FDWZ-five-trees" or params.method=="FDWZ-fixed-cubic" or params.method=="FDWZ-surface-first" or params.method=="FDWZ-normal-boundary" or params.method=="FDWZ-index-order" or params.method=="FDWZ-index-midpoint" or params.method=="QoZ-3D-GLL" or
            params.method.rfind("FDWZ-order-",0)==0 or params.method.rfind("FDWZ-reference-EB-",0)==0 or params.method.rfind("FDWZ-sampled-",0)==0 or params.method.rfind("FDWZ-position-bias-",0)==0)
            core::decompress(params);
        else
            baseline::decompress(params);
    }
}
