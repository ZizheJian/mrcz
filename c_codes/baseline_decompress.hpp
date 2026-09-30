#pragma once

#include "mrcz_core/mrcz.hpp"
#include "sz2_impl/sz2_decompress.tpp"
#include "sz3_impl/sz3_decompress.tpp"
#include "qoz_impl/qoz_decompress.tpp"
#include "zfp_impl/zfp_decompress.tpp"
#include "sperr_impl/sperr_decompress.tpp"
#include "mgard_impl/mgard_decompress.tpp"

namespace mrcz::baseline
{
    template<typename T>
    void decompress(param_class<T> &params)
    {
        auto backend=[]<typename U>(param_class<U> &backend_params,string field)
        {
            if (backend_params.method=="SZ2-1D" or backend_params.method=="SZ2-2D" or backend_params.method=="SZ2-3D-fused" or backend_params.method=="SZ2-3D-per-element" or backend_params.method=="SZ2-4D")
                sz2_decompress(backend_params,field);
            else if (backend_params.method=="SZ3-1D" or backend_params.method=="SZ3-2D" or backend_params.method=="SZ3-3D-fused" or backend_params.method=="SZ3-3D-per-element" or backend_params.method=="SZ3-4D")
                sz3_decompress(backend_params,field);
            else if (backend_params.method=="QoZ-1D" or backend_params.method=="QoZ-2D" or backend_params.method=="QoZ-3D-fused" or backend_params.method=="QoZ-3D-per-element" or backend_params.method=="QoZ-4D")
                qoz_decompress(backend_params,field);
            else if (backend_params.method=="ZFP-1D" or backend_params.method=="ZFP-2D" or backend_params.method=="ZFP-3D-fused" or backend_params.method=="ZFP-3D-per-element" or backend_params.method=="ZFP-4D")
                zfp_decompress(backend_params,field);
            else if (backend_params.method=="SPERR-2D" or backend_params.method=="SPERR-3D-fused" or backend_params.method=="SPERR-3D-per-element")
                sperr_decompress(backend_params,field);
            else if (backend_params.method.rfind("MGARD-",0)==0)
                mgard_decompress(backend_params,field);
            else
                throw runtime_error("未知的baseline解压方法: "+backend_params.method);
        };
        detail::decompress_file_with_backend(params,backend);
    }
}
