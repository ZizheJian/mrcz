#pragma once

#include "mrcz_core/mrcz.hpp"
#include "sz2_impl/sz2_1_compress.tpp"
#include "sz3_impl/sz3_1d_compress.tpp"
#include "sz3_impl/sz3_4d_compress.tpp"
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
            if (backend_params.method=="SZ2-3D-uniform" or backend_params.method=="SZ3-3D-uniform")
            {
                vector<vector<U>> backup_x=backend_params.x;
                for (int dim=0;dim<3;dim++)
                    if (backend_params.x[dim].size()>1)
                        for (uint32_t i=0;i<backend_params.x[dim].size();i++)
                            backend_params.x[dim][i]=(U)i/(backend_params.x[dim].size()-1)*2-1;
                    else
                        backend_params.x[dim][0]=0;
                if (backend_params.method=="SZ2-3D-uniform")
                    sz2_1_compress(backend_params,field);
                else
                    sz3_2_compress(backend_params,field);
                backend_params.x=backup_x;
            }
            else if (backend_params.method=="SZ2-3D-GLL")
                sz2_1_compress(backend_params,field);
            else if (backend_params.method=="SZ3-1D-uniform")
                sz3_1d_compress(backend_params,field);
            else if (backend_params.method=="SZ3-4D-uniform")
                sz3_4d_compress(backend_params,field);
            else if (backend_params.method=="QoZ-4D-uniform")
                qoz_compress(backend_params,field);
            else if (backend_params.method=="ZFP-4D-uniform")
                zfp_compress(backend_params,field);
            else if (backend_params.method=="SPERR-3D-uniform")
                sperr_compress(backend_params,field);
            else
                throw runtime_error("未知的baseline压缩方法: "+backend_params.method);
        };
        detail::compress_file_with_backend(params,backend);
    }
}
