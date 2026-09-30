#pragma once

#include "mrcz_core/mrcz.hpp"
#include "sz2_impl/sz2_compress.tpp"
#include "sz3_impl/sz3_compress.tpp"
#include "qoz_impl/qoz_compress.tpp"
#include "zfp_impl/zfp_compress.tpp"
#include "sperr_impl/sperr_compress.tpp"
#include "mgard_impl/mgard_compress.tpp"

namespace mrcz::baseline
{
    template<typename T>
    void compress(param_class<T> &params)
    {
        auto backend=[]<typename U>(param_class<U> &backend_params,string field)
        {
            if (backend_params.method=="SZ2-1D" or backend_params.method=="SZ2-2D" or backend_params.method=="SZ2-3D-fused" or backend_params.method=="SZ2-3D-per-element" or backend_params.method=="SZ2-4D")
                sz2_compress(backend_params,field);
            else if (backend_params.method=="SZ3-1D" or backend_params.method=="SZ3-2D" or backend_params.method=="SZ3-3D-fused" or backend_params.method=="SZ3-3D-per-element" or backend_params.method=="SZ3-4D")
                sz3_compress(backend_params,field);
            else if (backend_params.method=="QoZ-1D" or backend_params.method=="QoZ-2D" or backend_params.method=="QoZ-3D-fused" or backend_params.method=="QoZ-3D-per-element" or backend_params.method=="QoZ-4D")
                qoz_compress(backend_params,field);
            else if (backend_params.method=="ZFP-1D" or backend_params.method=="ZFP-2D" or backend_params.method=="ZFP-3D-fused" or backend_params.method=="ZFP-3D-per-element" or backend_params.method=="ZFP-4D")
                zfp_compress(backend_params,field);
            else if (backend_params.method=="SPERR-2D" or backend_params.method=="SPERR-3D-fused" or backend_params.method=="SPERR-3D-per-element")
                sperr_compress(backend_params,field);
            else if (backend_params.method.rfind("MGARD-",0)==0)
                mgard_compress(backend_params,field);
            else
                throw runtime_error("未知的baseline压缩方法: "+backend_params.method);
        };
        detail::compress_file_with_backend(params,backend);
    }
}
