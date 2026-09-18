#pragma once

#include <set>
#include "mrcz.hpp"

namespace mrcz::detail
{
    // 一维hierarchy记录每一层尚待细分的区间、最接近reference-space中点的预测目标，以及该层开始时已经重构的参考点。
    // 三个坐标方向分别生成hierarchy，之后再按照不同的维度顺序组合成三维预测plan。
    template<typename T>
    void fdwz_generate_hierarchy(vector<vector<uint32_t>> &left_indices,vector<vector<uint32_t>> &right_indices,vector<vector<uint32_t>> &mid_indices,
        vector<vector<uint32_t>> &compressed_indices,const vector<T> &x)
    {
        while (true)
        {
            mid_indices.push_back({});
            vector<uint32_t> new_left_indices,new_right_indices;
            for (uint32_t i=0;i<left_indices.back().size();i++)
            {
                uint32_t left_index=left_indices.back()[i],right_index=right_indices.back()[i];
                if(right_index-left_index<=1)
                    continue;
                uint32_t mid_index=left_index+1;
                T ideal_mid_pos=(x[left_index]+x[right_index])/2;
                // GLL点在reference space中不均匀，因此根据reference-space坐标选择中点，而不是直接使用索引中点。
                for (uint32_t j=left_index+2;j<right_index;j++)
                    if (fabs(x[j]-ideal_mid_pos)<fabs(x[mid_index]-ideal_mid_pos))
                        mid_index=j;
                mid_indices.back().push_back(mid_index);
                if (mid_index-left_index>1)
                    new_left_indices.push_back(left_index),new_right_indices.push_back(mid_index);
                if (right_index-mid_index>1)
                    new_left_indices.push_back(mid_index),new_right_indices.push_back(right_index);
            }
            if(mid_indices.back().empty())
            {
                mid_indices.pop_back();
                break;
            }
            set<uint32_t> temp_set;
            temp_set.insert(compressed_indices.back().begin(),compressed_indices.back().end());
            temp_set.insert(mid_indices.back().begin(),mid_indices.back().end());
            // 下一层只能引用角点和之前层级已经重构的点。
            compressed_indices.push_back(vector<uint32_t>(temp_set.begin(),temp_set.end()));
            if (new_left_indices.size()==0)
                break;
            else
                left_indices.push_back(new_left_indices),right_indices.push_back(new_right_indices);
        }
    }

    // 保存一个一维插值目标使用的最多四个参考点及其Lagrange权重。
    template<typename T>
    struct fdwz_interpolation_stencil
    {
        uint32_t point_num=0;
        uint32_t indices[4]={0,0,0,0};
        T weights[4]={0,0,0,0};
    };

    // 为一维hierarchy的每个目标生成linear或边界自适应的cubic插值stencil。
    template<typename T>
    vector<vector<fdwz_interpolation_stencil<T>>> fdwz_generate_interpolation_stencils(const vector<vector<uint32_t>> &left_indices,const vector<vector<uint32_t>> &right_indices,
        const vector<vector<uint32_t>> &mid_indices,const vector<vector<uint32_t>> &compressed_indices,const vector<T> &x,bool use_cubic)
    {
        // linear固定使用区间两端；cubic优先使用两侧各一个外点，边界处自动降为三点或两点插值。
        vector<vector<fdwz_interpolation_stencil<T>>> stencils(mid_indices.size());
        for (uint32_t level=0;level<mid_indices.size();level++)
        {
            stencils[level].resize(left_indices[level].size());
            for (uint32_t interval_id=0;interval_id<left_indices[level].size();interval_id++)
            {
                auto &stencil=stencils[level][interval_id];
                uint32_t left_index=left_indices[level][interval_id],right_index=right_indices[level][interval_id],mid_index=mid_indices[level][interval_id];
                stencil.point_num=2;
                stencil.indices[0]=left_index;
                stencil.indices[1]=right_index;
                uint32_t outer_left=left_index,outer_right=right_index;
                if (use_cubic)
                {
                    for (auto index:compressed_indices[level])
                    {
                        if (index<left_index and (outer_left==left_index or index>outer_left))
                            outer_left=index;
                        if (index>right_index and (outer_right==right_index or index<outer_right))
                            outer_right=index;
                    }
                    if (outer_left!=left_index and outer_right!=right_index)
                    {
                        stencil.point_num=4;
                        stencil.indices[0]=outer_left;
                        stencil.indices[1]=left_index;
                        stencil.indices[2]=right_index;
                        stencil.indices[3]=outer_right;
                    }
                    else if (outer_left!=left_index)
                    {
                        stencil.point_num=3;
                        stencil.indices[0]=outer_left;
                        stencil.indices[1]=left_index;
                        stencil.indices[2]=right_index;
                    }
                    else if (outer_right!=right_index)
                    {
                        stencil.point_num=3;
                        stencil.indices[0]=left_index;
                        stencil.indices[1]=right_index;
                        stencil.indices[2]=outer_right;
                    }
                }
                // 在reference-space GLL坐标上预计算Lagrange插值权重，压缩和解压必须使用同一组确定性权重。
                for (uint32_t i=0;i<stencil.point_num;i++)
                {
                    stencil.weights[i]=1;
                    for (uint32_t j=0;j<stencil.point_num;j++)
                        if (i!=j)
                            stencil.weights[i]*=(x[mid_index]-x[stencil.indices[j]])/(x[stencil.indices[i]]-x[stencil.indices[j]]);
                }
            }
        }
        return stencils;
    }

    // 返回一个element中所有不重复角点的C-order线性位置。
    inline vector<uint32_t> fdwz_generate_corner_positions(uint32_t n0,uint32_t n1,uint32_t n2)
    {
        // 长度为1的维度只有一个端点，因此1D/2D退化element不会生成重复角点。
        vector<uint32_t> indices0=n0==1?vector<uint32_t>{0}:vector<uint32_t>{0,n0-1};
        vector<uint32_t> indices1=n1==1?vector<uint32_t>{0}:vector<uint32_t>{0,n1-1};
        vector<uint32_t> indices2=n2==1?vector<uint32_t>{0}:vector<uint32_t>{0,n2-1};
        vector<uint32_t> positions;
        for (auto i0:indices0)
            for (auto i1:indices1)
                for (auto i2:indices2)
                    positions.push_back(i0*n1*n2+i1*n2+i2);
        return positions;
    }

    // 保存三维element内一个预测目标、参考点位置及对应权重。
    template<typename T>
    struct fdwz_prediction_stencil
    {
        // target和indices都是单个element内的C-order线性位置。
        uint32_t target=0;
        uint32_t point_num=0;
        uint32_t indices[4]={0,0,0,0};
        T weights[4]={0,0,0,0};
    };

    // 使用已经重构的element值和预计算stencil计算一个目标的预测值。
    template<typename T>
    inline T fdwz_predict_from_stencil(const T *element,const fdwz_prediction_stencil<T> &stencil)
    {
        T prediction=0;
        for(uint32_t i=0;i<stencil.point_num;i++)
            prediction+=stencil.weights[i]*element[stencil.indices[i]];
        return prediction;
    }

    // 按指定维度顺序把三个一维hierarchy展开为覆盖所有非角点的三维预测序列。
    template<typename T>
    vector<fdwz_prediction_stencil<T>> fdwz_generate_prediction_stencils(uint32_t n0,uint32_t n1,uint32_t n2,uint32_t n12,uint32_t max_level,
        const vector<vector<uint32_t>> &mid_indices0,const vector<vector<uint32_t>> &compressed_indices0,const vector<vector<fdwz_interpolation_stencil<T>>> &interpolation_stencils0,
        const vector<vector<uint32_t>> &mid_indices1,const vector<vector<uint32_t>> &compressed_indices1,const vector<vector<fdwz_interpolation_stencil<T>>> &interpolation_stencils1,
        const vector<vector<uint32_t>> &mid_indices2,const vector<vector<uint32_t>> &compressed_indices2,const vector<vector<fdwz_interpolation_stencil<T>>> &interpolation_stencils2,
        vector<uint32_t> &level_offsets,const vector<uint32_t> &dimension_order)
    {
        // 每个全局level依次推进三个维度。dimension_order决定同一level内哪个方向先产生可供后续方向引用的新点。
        vector<fdwz_prediction_stencil<T>> stencils;
        level_offsets={0};
        uint32_t level0=0,level1=0,level2=0;
        for (uint32_t level=0;level<max_level;level++)
        {
            for (auto dimension:dimension_order)
            {
                if (dimension==0 and n0>1 and level0<mid_indices0.size())
                {
                    for (uint32_t interval=0;interval<mid_indices0[level0].size();interval++)
                        for (auto i1:compressed_indices1[level1])
                            for (auto i2:compressed_indices2[level2])
                            {
                                uint32_t pos_base=i1*n2+i2;
                                fdwz_prediction_stencil<T> stencil;
                                stencil.target=pos_base+mid_indices0[level0][interval]*n12;
                                auto &interpolation_stencil=interpolation_stencils0[level0][interval];
                                stencil.point_num=interpolation_stencil.point_num;
                                for (uint32_t i=0;i<stencil.point_num;i++)
                                    stencil.indices[i]=pos_base+interpolation_stencil.indices[i]*n12,stencil.weights[i]=interpolation_stencil.weights[i];
                                stencils.push_back(stencil);
                            }
                    level0++;
                }
                else if (dimension==1 and n1>1 and level1<mid_indices1.size())
                {
                    for (auto i0:compressed_indices0[level0])
                        for (uint32_t interval=0;interval<mid_indices1[level1].size();interval++)
                            for (auto i2:compressed_indices2[level2])
                            {
                                uint32_t pos_base=i0*n12+i2;
                                fdwz_prediction_stencil<T> stencil;
                                stencil.target=pos_base+mid_indices1[level1][interval]*n2;
                                auto &interpolation_stencil=interpolation_stencils1[level1][interval];
                                stencil.point_num=interpolation_stencil.point_num;
                                for (uint32_t i=0;i<stencil.point_num;i++)
                                    stencil.indices[i]=pos_base+interpolation_stencil.indices[i]*n2,stencil.weights[i]=interpolation_stencil.weights[i];
                                stencils.push_back(stencil);
                            }
                    level1++;
                }
                else if (dimension==2 and n2>1 and level2<mid_indices2.size())
                {
                    for (auto i0:compressed_indices0[level0])
                        for (auto i1:compressed_indices1[level1])
                            for (uint32_t interval=0;interval<mid_indices2[level2].size();interval++)
                            {
                                uint32_t pos_base=i0*n12+i1*n2;
                                fdwz_prediction_stencil<T> stencil;
                                stencil.target=pos_base+mid_indices2[level2][interval];
                                auto &interpolation_stencil=interpolation_stencils2[level2][interval];
                                stencil.point_num=interpolation_stencil.point_num;
                                for (uint32_t i=0;i<stencil.point_num;i++)
                                    stencil.indices[i]=pos_base+interpolation_stencil.indices[i],stencil.weights[i]=interpolation_stencil.weights[i];
                                stencils.push_back(stencil);
                            }
                    level2++;
                }
            }
            level_offsets.push_back(stencils.size());
        }
        // plan必须恰好覆盖所有非角点，并且每个target只能出现一次；后面的residual_id依赖这个顺序。
        uint32_t n012=n0*n1*n2;
        vector<uint32_t> corner_positions=fdwz_generate_corner_positions(n0,n1,n2);
        if (stencils.size()!=n012-corner_positions.size())
            throw runtime_error("FDWZ prediction stencil数量错误");
        vector<bool> seen(n012,false);
        for (auto position:corner_positions)
            seen[position]=true;
        for (auto &stencil:stencils)
        {
            if (stencil.target>=n012 or seen[stencil.target])
                throw runtime_error("FDWZ prediction stencil target重复、越界或指向角点");
            seen[stencil.target]=true;
        }
        return stencils;
    }

    // 保存一个partition候选方案的完整预测序列和Huffman分层边界。
    template<typename T>
    struct fdwz_prediction_plan
    {
        // level_offsets既描述预测层边界，也用于把不同层的residual QB分配给多棵Huffman tree。
        vector<uint32_t> level_offsets;
        vector<fdwz_prediction_stencil<T>> stencils;
    };

    // 返回三维坐标方向的全部六种执行顺序。
    inline vector<vector<uint32_t>> fdwz_all_dimension_orders()
    {
        return {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    }

    // 集中保存由element shape和reference-space GLL坐标确定的FDWZ预测结构。
    template<typename T>
    struct fdwz_plan_set
    {
        vector<uint32_t> corner_positions;
        vector<fdwz_prediction_plan<T>> linear_plans,cubic_plans;
    };

    // 压缩端同时需要linear proxy和cubic plans；解压端只需要重建cubic plans。
    template<typename T>
    fdwz_plan_set<T> fdwz_build_plan_set(uint32_t n0,uint32_t n1,uint32_t n2,const vector<T> &x0,const vector<T> &x1,const vector<T> &x2,bool include_linear_plans)
    {
        uint32_t n12=n1*n2;
        vector<vector<uint32_t>> left_indices0{{0}},right_indices0{{n0-1}},mid_indices0,compressed_indices0;
        vector<vector<uint32_t>> left_indices1{{0}},right_indices1{{n1-1}},mid_indices1,compressed_indices1;
        vector<vector<uint32_t>> left_indices2{{0}},right_indices2{{n2-1}},mid_indices2,compressed_indices2;
        if(n0==1)
            compressed_indices0.push_back({0});
        else
        {
            compressed_indices0.push_back({0,n0-1});
            fdwz_generate_hierarchy(left_indices0,right_indices0,mid_indices0,compressed_indices0,x0);
        }
        if(n1==1)
            compressed_indices1.push_back({0});
        else
        {
            compressed_indices1.push_back({0,n1-1});
            fdwz_generate_hierarchy(left_indices1,right_indices1,mid_indices1,compressed_indices1,x1);
        }
        if(n2==1)
            compressed_indices2.push_back({0});
        else
        {
            compressed_indices2.push_back({0,n2-1});
            fdwz_generate_hierarchy(left_indices2,right_indices2,mid_indices2,compressed_indices2,x2);
        }
        uint32_t max_level=max(max(mid_indices0.size(),mid_indices1.size()),mid_indices2.size());
        fdwz_plan_set<T> plan_set;
        plan_set.corner_positions=fdwz_generate_corner_positions(n0,n1,n2);
        auto all_orders=fdwz_all_dimension_orders();
        auto build_plans=[&](bool use_cubic)
        {
            vector<fdwz_prediction_plan<T>> plans;
            auto interpolation_stencils0=fdwz_generate_interpolation_stencils(left_indices0,right_indices0,mid_indices0,compressed_indices0,x0,use_cubic);
            auto interpolation_stencils1=fdwz_generate_interpolation_stencils(left_indices1,right_indices1,mid_indices1,compressed_indices1,x1,use_cubic);
            auto interpolation_stencils2=fdwz_generate_interpolation_stencils(left_indices2,right_indices2,mid_indices2,compressed_indices2,x2,use_cubic);
            for(const auto &order:all_orders)
            {
                fdwz_prediction_plan<T> plan;
                plan.stencils=fdwz_generate_prediction_stencils(n0,n1,n2,n12,max_level,mid_indices0,compressed_indices0,interpolation_stencils0,
                    mid_indices1,compressed_indices1,interpolation_stencils1,mid_indices2,compressed_indices2,interpolation_stencils2,plan.level_offsets,order);
                if(!plans.empty() and plan.level_offsets!=plans[0].level_offsets)
                    throw runtime_error("FDWZ不同维度顺序的预测层结构不一致");
                plans.push_back(move(plan));
            }
            return plans;
        };
        if(include_linear_plans)
            plan_set.linear_plans=build_plans(false);
        plan_set.cubic_plans=build_plans(true);
        if(include_linear_plans and plan_set.linear_plans[0].level_offsets!=plan_set.cubic_plans[0].level_offsets)
            throw runtime_error("FDWZ linear和cubic plans的预测层结构不一致");
        return plan_set;
    }

}

