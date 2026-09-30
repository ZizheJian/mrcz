#pragma once

#include "baseline_compress.hpp"
#include "mrcz_core/mrcz.hpp"

namespace mrcz
{
    inline bool is_core_method(const string &method)
    {
        return method.empty() or method=="MRCZ" or method=="MRCZ-ablation" or method=="FDWZ" or method=="FDWZ-ablation";
    }

    inline bool is_baseline_method(const string &method)
    {
        return method=="SZ2-1D" or method=="SZ2-2D" or method=="SZ2-3D-fused" or method=="SZ2-3D-per-element" or method=="SZ2-4D" or
            method=="SZ3-1D" or method=="SZ3-2D" or method=="SZ3-3D-fused" or method=="SZ3-3D-per-element" or method=="SZ3-4D" or
            method=="QoZ-1D" or method=="QoZ-2D" or method=="QoZ-3D-fused" or method=="QoZ-3D-per-element" or method=="QoZ-4D" or
            method=="ZFP-1D" or method=="ZFP-2D" or method=="ZFP-3D-fused" or method=="ZFP-3D-per-element" or method=="ZFP-4D" or
            method=="SPERR-2D" or method=="SPERR-3D-fused" or method=="SPERR-3D-per-element" or
            method=="MGARD-1D" or method=="MGARD-2D" or method=="MGARD-3D-fused" or method=="MGARD-3D-per-element" or method=="MGARD-3D-per-element-GLL" or method=="MGARD-4D" or
            method=="MGARD-4D-GLL";
    }

    template<typename T>
    void compress(param_class<T> &params)
    {
        if (is_core_method(params.method))
            core::compress(params);
        else if (is_baseline_method(params.method))
            baseline::compress(params);
        else
            throw runtime_error("未知的压缩方法: "+params.method);
    }
}
