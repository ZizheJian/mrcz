#pragma once


#include "mrcz.hpp"

template<typename T>
struct sz3_normal_boundary_stencil
{
    uint32_t target=0;
    uint32_t reference_num=0;
    uint32_t references[8]={0,0,0,0,0,0,0,0};
    T weights[8]={0,0,0,0,0,0,0,0};
};

inline uint32_t sz3_normal_flatten(uint32_t *indices,uint32_t n1,uint32_t n2)
{
    return (indices[0]*n1+indices[1])*n2+indices[2];
}

template<typename T>
vector<sz3_normal_boundary_stencil<T>> sz3_generate_normal_boundary_stencils(vector<uint32_t> &shape,vector<vector<T>> &x)
{
    if (shape.size()!=3 or x.size()!=3)
        throw runtime_error("SZ3 normal-boundary要求三个空间维度");
    uint32_t n0=shape[0],n1=shape[1],n2=shape[2],n012=n0*n1*n2,active_dimension_num=0;
    for (auto n:shape)
        if (n>1)
        {
            active_dimension_num++;
            if (n<4)
                throw runtime_error("SZ3 normal-boundary要求每个有效维度至少有4个点");
        }
    if (active_dimension_num<2)
        throw runtime_error("SZ3 normal-boundary只支持二维或三维element");

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
                    known[sz3_normal_flatten(indices,n1,n2)]=true;
            }

    // XXX 偶数边长使用中间两个点，奇数边长使用中心及其两侧点，使两个法向都具有两个中心reference。
    vector<vector<uint32_t>> center_indices(3);
    for (uint32_t dim=0;dim<3;dim++)
    {
        if (shape[dim]==1)
            center_indices[dim]={0};
        else if (shape[dim]%2==0)
            center_indices[dim]={shape[dim]/2-1,shape[dim]/2};
        else
            center_indices[dim]={shape[dim]/2-1,shape[dim]/2,shape[dim]/2+1};
    }

    vector<sz3_normal_boundary_stencil<T>> stencils;
    // XXX 中心块只引用角点，使用真实GLL坐标上的双线性或三线性内插；中心块内部不互相引用。
    vector<sz3_normal_boundary_stencil<T>> center_stencils;
    for (auto i0:center_indices[0])
        for (auto i1:center_indices[1])
            for (auto i2:center_indices[2])
            {
                uint32_t target_indices[3]={i0,i1,i2};
                uint32_t target=sz3_normal_flatten(target_indices,n1,n2);
                if (known[target])
                    continue;
                sz3_normal_boundary_stencil<T> stencil;
                stencil.target=target;
                uint32_t corner_combination_num=1u<<active_dimension_num;
                for (uint32_t combination=0;combination<corner_combination_num;combination++)
                {
                    uint32_t corner_indices[3]={0,0,0},active_dimension=0;
                    T weight=1;
                    for (uint32_t dim=0;dim<3;dim++)
                    {
                        if (shape[dim]==1)
                            continue;
                        bool high=(combination>>active_dimension)&1u;
                        active_dimension++;
                        corner_indices[dim]=high?shape[dim]-1:0;
                        T low_weight=(x[dim].back()-x[dim][target_indices[dim]])/(x[dim].back()-x[dim].front());
                        weight*=high?1-low_weight:low_weight;
                    }
                    uint32_t reference=stencil.reference_num++;
                    stencil.references[reference]=sz3_normal_flatten(corner_indices,n1,n2);
                    stencil.weights[reference]=weight;
                }
                center_stencils.push_back(stencil);
            }
    for (auto &stencil:center_stencils)
        known[stencil.target]=true;
    stencils.insert(stencils.end(),center_stencils.begin(),center_stencils.end());

    auto append_available_boundary_rounds=[&]()
    {
        while (true)
        {
            vector<sz3_normal_boundary_stencil<T>> round_stencils;
            for (uint32_t i0=0;i0<n0;i0++)
                for (uint32_t i1=0;i1<n1;i1++)
                    for (uint32_t i2=0;i2<n2;i2++)
                    {
                        uint32_t target_indices[3]={i0,i1,i2};
                        uint32_t target=sz3_normal_flatten(target_indices,n1,n2),normal_num=0;
                        if (known[target])
                            continue;
                        for (uint32_t dim=0;dim<3;dim++)
                            if (shape[dim]>1 and (target_indices[dim]==0 or target_indices[dim]+1==shape[dim]))
                                normal_num++;
                        if (normal_num==0)
                            continue;

                        sz3_normal_boundary_stencil<T> stencil;
                        stencil.target=target;
                        bool all_normals_available=true;
                        for (uint32_t dim=0;dim<3;dim++)
                        {
                            if (shape[dim]==1 or (target_indices[dim]!=0 and target_indices[dim]+1!=shape[dim]))
                                continue;
                            uint32_t reference_coordinates[2]={0,0},found=0;
                            if (target_indices[dim]==0)
                            {
                                for (uint32_t coordinate=1;coordinate<shape[dim] and found<2;coordinate++)
                                {
                                    uint32_t reference_indices[3]={i0,i1,i2};
                                    reference_indices[dim]=coordinate;
                                    if (known[sz3_normal_flatten(reference_indices,n1,n2)])
                                        reference_coordinates[found++]=coordinate;
                                }
                            }
                            else
                            {
                                for (uint32_t coordinate=shape[dim]-1;coordinate>0 and found<2;)
                                {
                                    coordinate--;
                                    uint32_t reference_indices[3]={i0,i1,i2};
                                    reference_indices[dim]=coordinate;
                                    if (known[sz3_normal_flatten(reference_indices,n1,n2)])
                                        reference_coordinates[found++]=coordinate;
                                }
                            }
                            if (found<2)
                            {
                                all_normals_available=false;
                                break;
                            }
                            T target_x=x[dim][target_indices[dim]],x0=x[dim][reference_coordinates[0]],x1=x[dim][reference_coordinates[1]];
                            T weights[2]={(x1-target_x)/(x1-x0),(target_x-x0)/(x1-x0)};
                            for (uint32_t reference_id=0;reference_id<2;reference_id++)
                            {
                                uint32_t reference_indices[3]={i0,i1,i2};
                                reference_indices[dim]=reference_coordinates[reference_id];
                                uint32_t reference=stencil.reference_num++;
                                stencil.references[reference]=sz3_normal_flatten(reference_indices,n1,n2);
                                stencil.weights[reference]=weights[reference_id]/normal_num;
                            }
                        }
                        if (all_normals_available)
                            round_stencils.push_back(stencil);
                    }
            if (round_stencils.empty())
                break;
            for (auto &stencil:round_stencils)
                known[stencil.target]=true;
            stencils.insert(stencils.end(),round_stencils.begin(),round_stencils.end());
        }
    };
    append_available_boundary_rounds();

    // XXX 内部从中心向外扩展；轴向点先完成，之后面对角和体对角点平均多个方向的内插。
    uint32_t max_center_shell=max(max(n0,n1),n2);
    for (uint32_t center_shell=0;center_shell<2*max_center_shell;center_shell++)
        for (uint32_t outer_dimension_num=1;outer_dimension_num<=active_dimension_num;outer_dimension_num++)
        {
            vector<sz3_normal_boundary_stencil<T>> round_stencils;
            for (uint32_t i0=0;i0<n0;i0++)
                for (uint32_t i1=0;i1<n1;i1++)
                    for (uint32_t i2=0;i2<n2;i2++)
                    {
                        uint32_t target_indices[3]={i0,i1,i2};
                        bool is_boundary=false;
                        for (uint32_t dim=0;dim<3;dim++)
                            if (shape[dim]>1 and (target_indices[dim]==0 or target_indices[dim]+1==shape[dim]))
                                is_boundary=true;
                        if (is_boundary)
                            continue;
                        uint32_t target=sz3_normal_flatten(target_indices,n1,n2);
                        if (known[target])
                            continue;
                        uint32_t distances[3]={0,0,0},current_center_shell=0,current_outer_dimension_num=0;
                        for (uint32_t dim=0;dim<3;dim++)
                            if (shape[dim]>1)
                            {
                                distances[dim]=static_cast<uint32_t>(abs(static_cast<int64_t>(2*target_indices[dim])-static_cast<int64_t>(shape[dim]-1)));
                                current_center_shell=max(current_center_shell,distances[dim]);
                            }
                        for (uint32_t dim=0;dim<3;dim++)
                            if (shape[dim]>1 and distances[dim]==current_center_shell)
                                current_outer_dimension_num++;
                        if (current_center_shell!=center_shell or current_outer_dimension_num!=outer_dimension_num)
                            continue;

                        sz3_normal_boundary_stencil<T> stencil;
                        stencil.target=target;
                        bool all_directions_available=true;
                        for (uint32_t dim=0;dim<3;dim++)
                        {
                            if (shape[dim]==1 or distances[dim]!=center_shell)
                                continue;
                            uint32_t left_coordinate=target_indices[dim]-1,right_coordinate=target_indices[dim]+1;
                            uint32_t left_indices[3]={i0,i1,i2},right_indices[3]={i0,i1,i2};
                            left_indices[dim]=left_coordinate;
                            right_indices[dim]=right_coordinate;
                            uint32_t left=sz3_normal_flatten(left_indices,n1,n2),right=sz3_normal_flatten(right_indices,n1,n2);
                            while (!known[left] and left_coordinate>0)
                            {
                                left_coordinate--;
                                left_indices[dim]=left_coordinate;
                                left=sz3_normal_flatten(left_indices,n1,n2);
                            }
                            while (!known[right] and right_coordinate+1<shape[dim])
                            {
                                right_coordinate++;
                                right_indices[dim]=right_coordinate;
                                right=sz3_normal_flatten(right_indices,n1,n2);
                            }
                            if (!known[left] or !known[right])
                            {
                                all_directions_available=false;
                                break;
                            }
                            uint32_t left_reference=stencil.reference_num++,right_reference=stencil.reference_num++;
                            stencil.references[left_reference]=left;
                            stencil.references[right_reference]=right;
                            stencil.weights[left_reference]=(x[dim][right_coordinate]-x[dim][target_indices[dim]])/(x[dim][right_coordinate]-x[dim][left_coordinate])/outer_dimension_num;
                            stencil.weights[right_reference]=(x[dim][target_indices[dim]]-x[dim][left_coordinate])/(x[dim][right_coordinate]-x[dim][left_coordinate])/outer_dimension_num;
                        }
                        if (all_directions_available)
                            round_stencils.push_back(stencil);
                    }
            for (auto &stencil:round_stencils)
                known[stencil.target]=true;
            stencils.insert(stencils.end(),round_stencils.begin(),round_stencils.end());
            append_available_boundary_rounds();
        }

    uint32_t corner_num=1;
    for (auto n:shape)
        if (n>1)
            corner_num*=2;
    if (stencils.size()!=n012-corner_num)
        throw runtime_error("SZ3 normal-boundary预测表没有覆盖全部非角点");
    return stencils;
}

template<typename T>
inline T sz3_normal_boundary_predict(T *element,sz3_normal_boundary_stencil<T> &stencil)
{
    T prediction=0;
    for (uint32_t reference=0;reference<stencil.reference_num;reference++)
        prediction+=stencil.weights[reference]*element[stencil.references[reference]];
    return prediction;
}
