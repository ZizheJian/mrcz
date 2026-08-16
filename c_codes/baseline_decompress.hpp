#pragma once

#include "mrcz_core/mrcz.hpp"
#include "sz2_impl/sz2_1_decompress.tpp"
#include "zfp_impl/zfp_decompress.tpp"
#include "sperr_impl/sperr_decompress.tpp"

namespace mrcz::baseline
{
    template<typename T>
    void decompress(param_class<T> &params)
    {
        auto backend=[]<typename U>(param_class<U> &backend_params,const std::string &field)
        {
            if (backend_params.method=="sz2-0" or backend_params.method=="sz3-1")
            {
                std::vector<std::vector<U>> backup_x=backend_params.x;
                for (int dim=0;dim<3;dim++)
                    if (backend_params.x[dim].size()>1)
                        for (uint32_t i=0;i<backend_params.x[dim].size();i++)
                            backend_params.x[dim][i]=static_cast<U>(i)/(backend_params.x[dim].size()-1)*2-1;
                    else
                        backend_params.x[dim][0]=0;
                if (backend_params.method=="sz2-0")
                    sz2_1_decompress(backend_params,field);
                else
                    sz3_2_decompress(backend_params,field);
                backend_params.x=backup_x;
            }
            else if (backend_params.method=="sz2-1")
                sz2_1_decompress(backend_params,field);
            else if (backend_params.method=="zfp")
                zfp_decompress(backend_params,field);
            else if (backend_params.method=="sperr")
                sperr_decompress(backend_params,field);
            else
                throw std::runtime_error("未知的baseline解压方法: "+backend_params.method);
        };
        detail::decompress_file_with_backend(params,backend);
    }
}
