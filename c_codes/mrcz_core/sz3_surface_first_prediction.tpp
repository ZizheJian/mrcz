#pragma once


#include "mrcz.hpp"

template<typename T>
struct sz3_surface_first_stencil
{
    uint32_t target=0;
    uint32_t predictor_num=0;
    uint32_t left[3]={0,0,0};
    uint32_t right[3]={0,0,0};
    T left_weight[3]={0,0,0};
    T right_weight[3]={0,0,0};
};

inline uint32_t sz3_flatten_position(uint32_t i0,uint32_t i1,uint32_t i2,uint32_t n1,uint32_t n2)
{
    return (i0*n1+i1)*n2+i2;
}

template<typename T>
vector<sz3_surface_first_stencil<T>> sz3_generate_surface_first_stencils(vector<uint32_t> &shape,vector<vector<T>> &x)
{
    if (shape.size()!=3 or x.size()!=3)
        throw runtime_error("SZ3 surface-first要求三个空间维度");
    uint32_t n0=shape[0],n1=shape[1],n2=shape[2],n012=n0*n1*n2,active_dimension_num=0;
    for (auto n:shape)
        if (n>1)
            active_dimension_num++;
    if (active_dimension_num<2)
        throw runtime_error("SZ3 surface-first只支持二维或三维element");

    // XXX 角点已独立量化。之后按自由维数依次处理边、面、体内部。
    vector<bool> known(n012,false);
    for (uint32_t i0=0;i0<n0;i0++)
        for (uint32_t i1=0;i1<n1;i1++)
            for (uint32_t i2=0;i2<n2;i2++)
            {
                uint32_t indices[3]={i0,i1,i2};
                bool is_corner=true;
                for (uint32_t dim=0;dim<3;dim++)
                    if (shape[dim]>1 and indices[dim]!=0 and indices[dim]+1!=shape[dim])
                        is_corner=false;
                if (is_corner)
                    known[sz3_flatten_position(i0,i1,i2,n1,n2)]=true;
            }

    vector<sz3_surface_first_stencil<T>> stencils;
    uint32_t max_center_shell=max(max(n0,n1),n2);
    for (uint32_t free_dimension_num=1;free_dimension_num<=active_dimension_num;free_dimension_num++)
        for (uint32_t center_shell=0;center_shell<2*max_center_shell;center_shell++)
            // XXX 同一中心层先处理轴向点，再处理面对角和体对角点，避免reference依赖环。
            for (uint32_t outer_dimension_num=1;outer_dimension_num<=free_dimension_num;outer_dimension_num++)
            {
                vector<sz3_surface_first_stencil<T>> round_stencils;
                for (uint32_t i0=0;i0<n0;i0++)
                    for (uint32_t i1=0;i1<n1;i1++)
                        for (uint32_t i2=0;i2<n2;i2++)
                        {
                            uint32_t indices[3]={i0,i1,i2};
                            uint32_t distances[3]={0,0,0};
                            uint32_t current_free_dimension_num=0,current_center_shell=0,current_outer_dimension_num=0;
                            for (uint32_t dim=0;dim<3;dim++)
                                if (shape[dim]>1 and indices[dim]!=0 and indices[dim]+1!=shape[dim])
                                {
                                    current_free_dimension_num++;
                                    distances[dim]=static_cast<uint32_t>(abs(static_cast<int64_t>(2*indices[dim])-static_cast<int64_t>(shape[dim]-1)));
                                    current_center_shell=max(current_center_shell,distances[dim]);
                                }
                            for (uint32_t dim=0;dim<3;dim++)
                                if (shape[dim]>1 and indices[dim]!=0 and indices[dim]+1!=shape[dim] and distances[dim]==current_center_shell)
                                    current_outer_dimension_num++;
                            uint32_t target=sz3_flatten_position(i0,i1,i2,n1,n2);
                            if (known[target] or current_free_dimension_num!=free_dimension_num or current_center_shell!=center_shell or current_outer_dimension_num!=outer_dimension_num)
                                continue;

                            sz3_surface_first_stencil<T> stencil;
                            stencil.target=target;
                            for (uint32_t dim=0;dim<3;dim++)
                            {
                                if (shape[dim]==1 or indices[dim]==0 or indices[dim]+1==shape[dim] or distances[dim]!=center_shell)
                                    continue;
                                uint32_t left_coordinate=indices[dim]-1,right_coordinate=indices[dim]+1;
                                uint32_t left_indices[3]={i0,i1,i2},right_indices[3]={i0,i1,i2};
                                left_indices[dim]=left_coordinate;
                                right_indices[dim]=right_coordinate;
                                uint32_t left=sz3_flatten_position(left_indices[0],left_indices[1],left_indices[2],n1,n2);
                                uint32_t right=sz3_flatten_position(right_indices[0],right_indices[1],right_indices[2],n1,n2);
                                while (!known[left] and left_coordinate>0)
                                {
                                    left_coordinate--;
                                    left_indices[dim]=left_coordinate;
                                    left=sz3_flatten_position(left_indices[0],left_indices[1],left_indices[2],n1,n2);
                                }
                                while (!known[right] and right_coordinate+1<shape[dim])
                                {
                                    right_coordinate++;
                                    right_indices[dim]=right_coordinate;
                                    right=sz3_flatten_position(right_indices[0],right_indices[1],right_indices[2],n1,n2);
                                }
                                if (!known[left] or !known[right])
                                    throw runtime_error("SZ3 surface-first预测表存在未建立的reference");
                                uint32_t predictor=stencil.predictor_num++;
                                stencil.left[predictor]=left;
                                stencil.right[predictor]=right;
                                stencil.left_weight[predictor]=(x[dim][right_coordinate]-x[dim][indices[dim]])/(x[dim][right_coordinate]-x[dim][left_coordinate]);
                                stencil.right_weight[predictor]=(x[dim][indices[dim]]-x[dim][left_coordinate])/(x[dim][right_coordinate]-x[dim][left_coordinate]);
                            }
                            if (stencil.predictor_num!=outer_dimension_num)
                                throw runtime_error("SZ3 surface-first预测方向数量错误");
                            round_stencils.push_back(stencil);
                        }
                // XXX 同一子轮内不互相引用，保证偶数边长的中心点和对称位置使用对称reference。
                for (auto &stencil:round_stencils)
                    known[stencil.target]=true;
                stencils.insert(stencils.end(),round_stencils.begin(),round_stencils.end());
            }
    uint32_t corner_num=1;
    for (auto n:shape)
        if (n>1)
            corner_num*=2;
    if (stencils.size()!=n012-corner_num)
        throw runtime_error("SZ3 surface-first预测表没有覆盖全部非角点");
    return stencils;
}

template<typename T>
inline T sz3_surface_first_predict(T *element,sz3_surface_first_stencil<T> &stencil)
{
    T prediction=0;
    for (uint32_t predictor=0;predictor<stencil.predictor_num;predictor++)
        prediction+=stencil.left_weight[predictor]*element[stencil.left[predictor]]+stencil.right_weight[predictor]*element[stencil.right[predictor]];
    return prediction/stencil.predictor_num;
}
