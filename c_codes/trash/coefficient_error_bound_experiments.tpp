#pragma once

// XXX 已淘汰的MRCZ分模式系数EB实验，不参与正式版本编译。
// XXX L1、RMS、product-0.5和总阶数方案在三个小数据上均未获得稳定的同PSNR BR收益，结果见TODO.md的历史记录。

namespace mrcz::detail
{
    inline bool is_order_coefficient_error_bound_method(string method)
    {
        return method=="MRCZ-3D-GLL-coeff-order-exponential-coarse" or method=="MRCZ-3D-GLL-coeff-order-exponential-fine" or method=="MRCZ-3D-GLL-coeff-order-polynomial-coarse" or method=="MRCZ-3D-GLL-coeff-order-polynomial-fine";
    }

    template<typename T>
    vector<T> compute_mode_coefficient_error_bounds(vector<vector<T>> &neg_L,vector<uint32_t> &element_shape,T uniform_eb,string norm)
    {
        vector<vector<T>> one_dimensional_sensitivity(3);
        for (int dim=0;dim<3;dim++)
        {
            uint32_t n=element_shape[dim];
            one_dimensional_sensitivity[dim].resize(n);
            for (uint32_t mode=0;mode<n;mode++)
            {
                double sum=0;
                for (uint32_t node=0;node<n;node++)
                {
                    double value=neg_L[dim][mode*n+node];
                    sum+=norm=="L1" ? fabs(value) : value*value;
                }
                one_dimensional_sensitivity[dim][mode]=norm=="L1" ? sum/n : sqrt(sum/n);
            }
        }

        uint32_t n0=element_shape[0],n1=element_shape[1],n2=element_shape[2],n12=n1*n2,n012=n0*n12;
        vector<T> sensitivity(n012),coefficient_eb(n012);
        double normalization_sum=0;
        for (uint32_t i=0;i<n0;i++)
            for (uint32_t j=0;j<n1;j++)
                for (uint32_t k=0;k<n2;k++)
                {
                    uint32_t index=i*n12+j*n2+k;
                    sensitivity[index]=one_dimensional_sensitivity[0][i]*one_dimensional_sensitivity[1][j]*one_dimensional_sensitivity[2][k];
                    normalization_sum+=norm=="L1" ? sensitivity[index] : sensitivity[index]*sensitivity[index];
                }

        T normalization=norm=="L1" ? normalization_sum/n012 : sqrt(normalization_sum/n012);
        for (uint32_t index=0;index<n012;index++)
        {
            if (sensitivity[index]<=0)
                throw runtime_error("反变换模式的敏感度必须大于0");
            coefficient_eb[index]=uniform_eb*normalization/sensitivity[index];
        }
        return coefficient_eb;
    }

    template<typename T>
    vector<T> compute_product_coefficient_error_bounds(vector<uint32_t> &element_shape,T uniform_eb,bool high_order_coarse)
    {
        uint32_t n0=element_shape[0],n1=element_shape[1],n2=element_shape[2],n12=n1*n2,n012=n0*n12;
        vector<T> coefficient_eb(n012);
        for (uint32_t i=0;i<n0;i++)
            for (uint32_t j=0;j<n1;j++)
                for (uint32_t k=0;k<n2;k++)
                {
                    uint32_t index=i*n12+j*n2+k;
                    T factor=sqrt(static_cast<T>(i+1)*static_cast<T>(j+1)*static_cast<T>(k+1));
                    coefficient_eb[index]=high_order_coarse ? uniform_eb*factor : uniform_eb/factor;
                }
        return coefficient_eb;
    }

    template<typename T>
    vector<T> compute_order_coefficient_error_bounds(vector<uint32_t> &element_shape,T uniform_eb,bool exponential,bool high_order_coarse)
    {
        uint32_t n0=element_shape[0],n1=element_shape[1],n2=element_shape[2],n12=n1*n2,n012=n0*n12;
        uint32_t low_order_max=1,max_order=n0+n1+n2-3;
        vector<T> coefficient_eb(n012,uniform_eb);
        if (max_order<=low_order_max)
            return coefficient_eb;
        for (uint32_t i=0;i<n0;i++)
            for (uint32_t j=0;j<n1;j++)
                for (uint32_t k=0;k<n2;k++)
                {
                    uint32_t order=i+j+k,index=i*n12+j*n2+k;
                    if (order<=low_order_max)
                        continue;
                    T normalized_order=static_cast<T>(order-low_order_max)/(max_order-low_order_max);
                    T factor=exponential ? pow(T(4),normalized_order) : 1+3*normalized_order*normalized_order;
                    coefficient_eb[index]=high_order_coarse ? uniform_eb*factor : uniform_eb/factor;
                }
        return coefficient_eb;
    }
}
