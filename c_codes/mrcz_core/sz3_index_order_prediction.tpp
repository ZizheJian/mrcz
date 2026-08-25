#pragma once


#include "mrcz.hpp"

template<typename T>
struct sz3_index_order_stencil
{
    uint32_t target=0;
    uint32_t reference0=0;
    uint32_t reference1=0;
    T weight0=1;
    T weight1=0;
    uint32_t reference_num=1;
};

template<typename T>
vector<sz3_index_order_stencil<T>> sz3_generate_index_order_stencils(vector<uint32_t> &shape,vector<vector<T>> &x)
{
    if (shape.size()!=3 or x.size()!=3)
        throw runtime_error("SZ3 index-order要求三维element和三组坐标");
    uint32_t n0=shape[0],n1=shape[1],n2=shape[2],n12=n1*n2,n012=n0*n1*n2;
    uint32_t max_size=max(n0,max(n1,n2)),max_level=0;
    while ((uint64_t(1)<<max_level)<max_size)
        max_level++;
    vector<sz3_index_order_stencil<T>> stencils;
    vector<bool> known(n012,false);
    known[0]=true;
    auto flat_index=[&](uint32_t i0,uint32_t i1,uint32_t i2){return i0*n12+i1*n2+i2;};
    auto append_stencil=[&](uint32_t i0,uint32_t i1,uint32_t i2,uint32_t dimension,uint32_t stride,uint32_t line_size)
    {
        uint32_t target_indices[3]={i0,i1,i2};
        uint32_t reference0_indices[3]={i0,i1,i2};
        uint32_t reference1_indices[3]={i0,i1,i2};
        uint32_t line_index=target_indices[dimension]/stride;
        sz3_index_order_stencil<T> stencil;
        stencil.target=flat_index(target_indices[0],target_indices[1],target_indices[2]);
        if (line_index+1<line_size)
        {
            reference0_indices[dimension]=(line_index-1)*stride;
            reference1_indices[dimension]=(line_index+1)*stride;
            stencil.reference_num=2;
        }
        else if (line_size<4)
            reference0_indices[dimension]=(line_index-1)*stride;
        else
        {
            reference0_indices[dimension]=(line_index-3)*stride;
            reference1_indices[dimension]=(line_index-1)*stride;
            stencil.reference_num=2;
        }
        stencil.reference0=flat_index(reference0_indices[0],reference0_indices[1],reference0_indices[2]);
        if (stencil.reference_num==2)
        {
            stencil.reference1=flat_index(reference1_indices[0],reference1_indices[1],reference1_indices[2]);
            T target_x=x[dimension][target_indices[dimension]],reference0_x=x[dimension][reference0_indices[dimension]],reference1_x=x[dimension][reference1_indices[dimension]];
            stencil.weight0=(reference1_x-target_x)/(reference1_x-reference0_x);
            stencil.weight1=(target_x-reference0_x)/(reference1_x-reference0_x);
        }
        if (stencil.target>=n012 or known[stencil.target] or !known[stencil.reference0] or (stencil.reference_num==2 and !known[stencil.reference1]))
            throw runtime_error("SZ3 index-order生成了重复target或尚未重建的reference");
        known[stencil.target]=true;
        stencils.push_back(stencil);
    };
    for (uint32_t level=max_level;level>0;level--)
    {
        uint32_t stride=uint32_t(1)<<(level-1),double_stride=2*stride;
        for (uint32_t dimension=0;dimension<3;dimension++)
        {
            uint32_t line_size=(shape[dimension]-1)/stride+1;
            if (line_size<=1)
                continue;
            uint32_t begin[3]={0,0,0},end[3]={n0,n1,n2},step[3]={dimension>0?stride:double_stride,dimension>1?stride:double_stride,double_stride};
            begin[dimension]=stride;
            end[dimension]=(line_size-1)*stride;
            step[dimension]=double_stride;
            for (uint32_t i0=begin[0];i0<end[0];i0+=step[0])
                for (uint32_t i1=begin[1];i1<end[1];i1+=step[1])
                    for (uint32_t i2=begin[2];i2<end[2];i2+=step[2])
                        append_stencil(i0,i1,i2,dimension,stride,line_size);
            if (line_size%2==0)
            {
                begin[dimension]=(line_size-1)*stride;
                end[dimension]=begin[dimension]+1;
                step[dimension]=1;
                for (uint32_t i0=begin[0];i0<end[0];i0+=step[0])
                    for (uint32_t i1=begin[1];i1<end[1];i1+=step[1])
                        for (uint32_t i2=begin[2];i2<end[2];i2+=step[2])
                            append_stencil(i0,i1,i2,dimension,stride,line_size);
            }
        }
    }
    if (stencils.size()!=n012-1)
        throw runtime_error("SZ3 index-order没有覆盖除原点外的所有节点");
    return stencils;
}

template<typename T>
inline T sz3_index_order_predict(T *element,sz3_index_order_stencil<T> &stencil)
{
    T prediction=stencil.weight0*element[stencil.reference0];
    if (stencil.reference_num==2)
        prediction+=stencil.weight1*element[stencil.reference1];
    return prediction;
}
