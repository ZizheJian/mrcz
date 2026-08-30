#pragma once

#include "mrcz_core/mrcz.hpp"
#include "sz2_impl/sz2_api_compress.tpp"
#include "sz3_impl/sz3_api_compress.tpp"
#include "qoz_impl/qoz_compress.tpp"
#include "zfp_impl/zfp_compress.tpp"
#include "sperr_impl/sperr_compress.tpp"

namespace mrcz::baseline
{
    template<typename T>
    void compress(param_class<T> &params)
    {
        auto backend=[]<typename U>(param_class<U> &backend_params,string field)
        {
            if (backend_params.method=="SZ2-1D-uniform" or backend_params.method=="SZ2-2D-uniform" or backend_params.method=="SZ2-3D-uniform-fused" or backend_params.method=="SZ2-3D-uniform-per-element" or backend_params.method=="SZ2-4D-uniform")
                sz2_api_compress(backend_params,field);
            else if (backend_params.method=="SZ3-1D-uniform" or backend_params.method=="SZ3-2D-uniform" or backend_params.method=="SZ3-3D-uniform-fused" or backend_params.method=="SZ3-3D-uniform-per-element" or backend_params.method=="SZ3-4D-uniform")
                sz3_api_compress(backend_params,field);
            else if (backend_params.method=="QoZ-1D-uniform" or backend_params.method=="QoZ-2D-uniform" or backend_params.method=="QoZ-3D-uniform-fused" or backend_params.method=="QoZ-3D-uniform-per-element" or backend_params.method=="QoZ-4D-uniform")
                qoz_compress(backend_params,field);
            else if (backend_params.method=="ZFP-1D-uniform" or backend_params.method=="ZFP-2D-uniform" or backend_params.method=="ZFP-3D-uniform-fused" or backend_params.method=="ZFP-3D-uniform-per-element" or backend_params.method=="ZFP-4D-uniform")
                zfp_compress(backend_params,field);
            else if (backend_params.method=="SPERR-2D-uniform" or backend_params.method=="SPERR-3D-uniform-fused" or backend_params.method=="SPERR-3D-uniform-per-element")
                sperr_compress(backend_params,field);
            else
                throw runtime_error("未知的baseline压缩方法: "+backend_params.method);
        };
        detail::compress_file_with_backend(params,backend);
    }
}
