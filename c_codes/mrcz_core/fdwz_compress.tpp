#pragma once

#include <array>
#include <memory>
#include <set>
#include "mrcz.hpp"
#include "huffman.hpp"
#include "quantize.tpp"
#include <zstd.h>

namespace mrcz::detail
{
    // 一维hierarchy记录每一层尚待细分的区间、最接近物理中点的预测目标，以及该层开始时已经重构的参考点。
    // 三个坐标方向分别生成hierarchy，之后再按照不同的维度顺序组合成三维预测plan。
    template<typename T>
    void fdwz_generate_hierarchy(vector<vector<uint32_t>> &left_indices,vector<vector<uint32_t>> &right_indices,vector<vector<uint32_t>> &mid_indices,
        vector<vector<uint32_t>> &compressed_indices,vector<T> &x)
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
                // GLL点不均匀，因此使用物理坐标选择中点，而不是直接使用索引中点。
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
    vector<vector<fdwz_interpolation_stencil<T>>> fdwz_generate_interpolation_stencils(vector<vector<uint32_t>> &left_indices,vector<vector<uint32_t>> &right_indices,
        vector<vector<uint32_t>> &mid_indices,vector<vector<uint32_t>> &compressed_indices,vector<T> &x,bool use_cubic)
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
                // 在真实GLL坐标上预计算Lagrange插值权重，压缩和解压必须使用同一组确定性权重。
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
    inline T fdwz_predict_from_stencil(T *element,fdwz_prediction_stencil<T> &stencil)
    {
        T prediction=0;
        for(uint32_t i=0;i<stencil.point_num;i++)
            prediction+=stencil.weights[i]*element[stencil.indices[i]];
        return prediction;
    }

    // 按指定维度顺序把三个一维hierarchy展开为覆盖所有非角点的三维预测序列。
    template<typename T>
    vector<fdwz_prediction_stencil<T>> fdwz_generate_prediction_stencils(uint32_t n0,uint32_t n1,uint32_t n2,uint32_t n12,uint32_t max_level,
        vector<vector<uint32_t>> &mid_indices0,vector<vector<uint32_t>> &compressed_indices0,vector<vector<fdwz_interpolation_stencil<T>>> &interpolation_stencils0,
        vector<vector<uint32_t>> &mid_indices1,vector<vector<uint32_t>> &compressed_indices1,vector<vector<fdwz_interpolation_stencil<T>>> &interpolation_stencils1,
        vector<vector<uint32_t>> &mid_indices2,vector<vector<uint32_t>> &compressed_indices2,vector<vector<fdwz_interpolation_stencil<T>>> &interpolation_stencils2,
        vector<uint32_t> &level_offsets,vector<uint32_t> &dimension_order)
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

    // 复用候选方案评估所需的临时数组，避免每个partition连续评估6个方案时反复申请大块内存。
    template<typename T>
    struct fdwz_sampling_workspace
    {
        vector<T> original_sample,sample;
        vector<vector<int32_t>> tree_qb;
        vector<unsigned char> huffman_bs,zstd_bs;
    };

    // INT32_MIN作为unpredictable标记；超出QB范围或浮点重构无法严格满足误差界时原值无损保存。
    template<bool check_qb_range,bool strict_error_check,typename T>
    inline int32_t fdwz_quantize_or_save_unpredictable(T &value,T prediction,T abs_eb,T neg_2_eb,vector<T> &unpredictable)
    {
        T scaled=(prediction-value)*neg_2_eb;
        T rounded=round(scaled);
        if constexpr(!check_qb_range and !strict_error_check)
        {
            int32_t qb=static_cast<int32_t>(rounded);
            value=prediction-qb*2*abs_eb;
            return qb;
        }
        else
        {
            long double rounded_bound=static_cast<long double>(rounded);
            if(!check_qb_range or (rounded_bound>=-static_cast<long double>(numeric_limits<int32_t>::max()) and rounded_bound<=static_cast<long double>(numeric_limits<int32_t>::max())))
            {
                int32_t qb=static_cast<int32_t>(rounded);
                T reconstructed=prediction-qb*2*abs_eb;
                if constexpr(!strict_error_check)
                {
                    value=reconstructed;
                    return qb;
                }
                else if(isfinite(reconstructed) and abs(reconstructed-value)<=abs_eb)
                {
                    value=reconstructed;
                    return qb;
                }
            }
        }
        unpredictable.push_back(value);
        return numeric_limits<int32_t>::min();
    }

    // 相对EB和stencil权重可以在field级证明全部残差QB安全；角点仍逐个检查其相对零的绝对QB。
    template<typename T>
    bool fdwz_residual_qb_range_is_safe(T data_min,T data_max,T abs_eb,vector<fdwz_prediction_plan<T>> &plans)
    {
        long double lower=data_min,upper=data_max,range=upper-lower,max_abs=max(abs(lower),abs(upper));
        long double epsilon=numeric_limits<T>::epsilon();
        for(auto &plan:plans)
            for(auto &stencil:plan.stencils)
            {
                long double weight_sum=0,weight_abs_sum=0;
                for(uint32_t i=0;i<stencil.point_num;i++)
                    weight_sum+=stencil.weights[i],weight_abs_sum+=abs(static_cast<long double>(stencil.weights[i]));
                long double residual_bound=(range+2*abs_eb)*(weight_abs_sum+1)+max_abs*(abs(weight_sum-1)+64*epsilon*weight_abs_sum);
                long double qb_bound=residual_bound/(2*static_cast<long double>(abs_eb))+0.5L;
                if(!isfinite(qb_bound) or qb_bound>numeric_limits<int32_t>::max())
                    return false;
            }
        return true;
    }

    // 返回三维坐标方向的全部六种执行顺序。
    inline vector<vector<uint32_t>> fdwz_all_dimension_orders()
    {
        return {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    }

    // 在抽样element上执行预测、量化、Huffman和Zstd，以最终字节数评价一个候选plan。
    template<typename T>
    uint64_t fdwz_sample_compressed_size(vector<T> &original_sample,uint32_t n012,vector<uint32_t> &corner_positions,fdwz_prediction_plan<T> &plan,
        vector<uint32_t> &residual_group_offsets,uint32_t huffman_tree_num,T abs_eb,unsigned char configuration,fdwz_sampling_workspace<T> &workspace)
    {
        uint64_t sample_num=original_sample.size()/n012;
        uint32_t corner_num=corner_positions.size();
        uint32_t residual_num=n012-corner_num;
        T neg_2_eb=1/(2*abs_eb);
        // 预测会原地写入重构值；每个候选plan必须从完全相同的原始sample开始。
        auto &sample=workspace.sample;
        sample=original_sample;
        // 直接量化到Huffman tree需要的最终布局，避免每个候选方案再完整重排一遍QB。
        auto &tree_qb=workspace.tree_qb;
        vector<T> unpredictable;
        tree_qb.resize(huffman_tree_num);
        if (huffman_tree_num==1)
            tree_qb[0].resize(sample_num*n012);
        else
        {
            tree_qb[0].resize(sample_num*corner_num);
            for (uint32_t group_id=0;group_id+1<residual_group_offsets.size();group_id++)
                tree_qb[group_id+1].resize(sample_num*(residual_group_offsets[group_id+1]-residual_group_offsets[group_id]));
        }
        for (uint64_t sample_id=0;sample_id<sample_num;sample_id++)
        {
            T *element=sample.data()+sample_id*n012;
            int32_t *element_corner_qb=tree_qb[0].data()+sample_id*(huffman_tree_num==1?n012:corner_num);
            for (uint32_t corner_id=0;corner_id<corner_num;corner_id++)
            {
                uint32_t target=corner_positions[corner_id];
                element_corner_qb[corner_id]=fdwz_quantize_or_save_unpredictable<true,false>(element[target],T(0),abs_eb,neg_2_eb,unpredictable);
            }
            if (huffman_tree_num==1)
            {
                int32_t *element_residual_qb=element_corner_qb+corner_num;
                for (uint32_t residual_id=0;residual_id<residual_num;residual_id++)
                {
                    auto &stencil=plan.stencils[residual_id];
                    T prediction=fdwz_predict_from_stencil(element,stencil);
                    element_residual_qb[residual_id]=fdwz_quantize_or_save_unpredictable<true,false>(element[stencil.target],prediction,abs_eb,neg_2_eb,unpredictable);
                }
            }
            else
            {
                for (uint32_t group_id=0;group_id+1<residual_group_offsets.size();group_id++)
                {
                    uint32_t group_begin=residual_group_offsets[group_id],group_size=residual_group_offsets[group_id+1]-group_begin;
                    int32_t *element_group_qb=tree_qb[group_id+1].data()+sample_id*group_size;
                    for (uint32_t local_id=0;local_id<group_size;local_id++)
                    {
                        auto &stencil=plan.stencils[group_begin+local_id];
                        T prediction=fdwz_predict_from_stencil(element,stencil);
                        element_group_qb[local_id]=fdwz_quantize_or_save_unpredictable<true,false>(element[stencil.target],prediction,abs_eb,neg_2_eb,unpredictable);
                    }
                }
            }
        }
        // sample评分使用保守的临时空间；首字节保存当前候选plan id，供Zstd评分包含相同的metadata开销。
        auto &huffman_bs_vec=workspace.huffman_bs;
        huffman_bs_vec.resize(sample_num*n012*32+4096*huffman_tree_num+sizeof(uint64_t)+1);
        unsigned char *huffman_bs=huffman_bs_vec.data();
        *huffman_bs++=configuration|0xc0;
        for (auto &values:tree_qb)
        {
            mrcz::huffman_encoder_class<int32_t> encoder;
            encoder.preprocess_encode(values.data(),values.size());
            encoder.save(huffman_bs);
            encoder.encode(values.data(),values.size(),huffman_bs);
        }
        if(!unpredictable.empty())
        {
            memcpy(huffman_bs,unpredictable.data(),unpredictable.size()*sizeof(T));
            huffman_bs+=unpredictable.size()*sizeof(T);
        }
        uint64_t huffman_bs_len=huffman_bs-huffman_bs_vec.data();
        auto &zstd_bs=workspace.zstd_bs;
        zstd_bs.resize(ZSTD_compressBound(huffman_bs_len));
        uint64_t zstd_bs_len=ZSTD_compress(zstd_bs.data(),zstd_bs.size(),huffman_bs_vec.data(),huffman_bs_len,3);
        if (ZSTD_isError(zstd_bs_len))
            throw runtime_error("ZSTD_compress error while sampling FDWZ configuration");
        return zstd_bs_len;
    }

    // 压缩一个partition；QB固定为int32，plan字节的bit 7/6分别标记32位QB和unpredictable原值流。
    template<bool measure_phase_time,bool check_residual_qb_range,bool strict_error_check,typename T>
    void fdwz_compress_block(mrcz::param_class<T> &params,uint64_t block_id,uint32_t n012,uint32_t corner_num,uint32_t residual_num,vector<uint32_t> &corner_positions,
        vector<fdwz_prediction_stencil<T>> &prediction_stencils,vector<uint32_t> &residual_group_offsets,array<uint32_t,5> &tree_value_num,uint32_t huffman_tree_num,T abs_eb,T neg_2_eb,
        unsigned char selected_plan_id,int32_t *block_qb,vector<uint64_t> &unpredictable_counts,vector<double> &prediction_quantization_time,vector<double> &huffman_time,vector<double> &zstd_time)
    {
        auto &block=params.partitions[block_id];
        auto &data=params.data;
        uint64_t begin_index=block.begin_index,end_index=block.end_index,element_num=end_index-begin_index;
        array<uint64_t,6> tree_offsets{};
        for(uint32_t tree_id=0;tree_id<huffman_tree_num;tree_id++)
            tree_offsets[tree_id+1]=tree_offsets[tree_id]+tree_value_num[tree_id]*element_num;
        double phase_start_time=0;
        if constexpr(measure_phase_time)
            phase_start_time=omp_get_wtime();
        vector<T> unpredictable;
        for(uint64_t i=begin_index;i<end_index;i++)
        {
            uint64_t relative_element_id=i-begin_index;
            int32_t *element_corner_qb=block_qb+relative_element_id*(huffman_tree_num==1?n012:corner_num);
            for(uint32_t corner_id=0;corner_id<corner_num;corner_id++)
            {
                uint32_t target=corner_positions[corner_id];
                element_corner_qb[corner_id]=fdwz_quantize_or_save_unpredictable<true,strict_error_check>(data[i*n012+target],T(0),abs_eb,neg_2_eb,unpredictable);
            }
        }
        if(huffman_tree_num==1)
        {
            for(uint64_t i=begin_index;i<end_index;i++)
            {
                T *element=data.data()+i*n012;
                uint64_t relative_element_id=i-begin_index;
                int32_t *element_residual_qb=block_qb+relative_element_id*n012+corner_num;
                for(uint32_t residual_id=0;residual_id<residual_num;residual_id++)
                {
                    auto &stencil=prediction_stencils[residual_id];
                    T prediction=fdwz_predict_from_stencil(element,stencil);
                    element_residual_qb[residual_id]=fdwz_quantize_or_save_unpredictable<check_residual_qb_range,strict_error_check>(element[stencil.target],prediction,abs_eb,neg_2_eb,unpredictable);
                }
            }
        }
        else
        {
            constexpr uint64_t batch_size=4;
            uint64_t batch_end=begin_index+(element_num/batch_size)*batch_size;
            for(uint64_t batch_begin=begin_index;batch_begin<batch_end;batch_begin+=batch_size)
            {
                uint64_t batch_relative_element_id=batch_begin-begin_index;
                for(uint32_t group_id=0;group_id+1<residual_group_offsets.size();group_id++)
                {
                    uint32_t group_begin=residual_group_offsets[group_id],group_size=residual_group_offsets[group_id+1]-group_begin;
                    int32_t *batch_group_qb=block_qb+tree_offsets[group_id+1]+batch_relative_element_id*group_size;
                    for(uint32_t local_id=0;local_id<group_size;local_id++)
                    {
                        auto &stencil=prediction_stencils[group_begin+local_id];
                        #pragma GCC unroll 4
                        for(uint64_t lane=0;lane<batch_size;lane++)
                        {
                            T *batch_element=data.data()+(batch_begin+lane)*n012;
                            T prediction=fdwz_predict_from_stencil(batch_element,stencil);
                            batch_group_qb[lane*group_size+local_id]=fdwz_quantize_or_save_unpredictable<check_residual_qb_range,strict_error_check>(batch_element[stencil.target],prediction,abs_eb,neg_2_eb,unpredictable);
                        }
                    }
                }
            }
            for(uint64_t i=batch_end;i<end_index;i++)
            {
                T *element=data.data()+i*n012;
                uint64_t relative_element_id=i-begin_index;
                for(uint32_t group_id=0;group_id+1<residual_group_offsets.size();group_id++)
                {
                    uint32_t group_begin=residual_group_offsets[group_id],group_size=residual_group_offsets[group_id+1]-group_begin;
                    int32_t *element_group_qb=block_qb+tree_offsets[group_id+1]+relative_element_id*group_size;
                    for(uint32_t local_id=0;local_id<group_size;local_id++)
                    {
                        auto &stencil=prediction_stencils[group_begin+local_id];
                        T prediction=fdwz_predict_from_stencil(element,stencil);
                        element_group_qb[local_id]=fdwz_quantize_or_save_unpredictable<check_residual_qb_range,strict_error_check>(element[stencil.target],prediction,abs_eb,neg_2_eb,unpredictable);
                    }
                }
            }
        }
        if constexpr(measure_phase_time)
            prediction_quantization_time[block_id]=omp_get_wtime()-phase_start_time;
        unpredictable_counts[block_id]=unpredictable.size();
        uint64_t huffman_bs_capacity=n012*element_num*sizeof(int64_t)+unpredictable.size()*sizeof(T)+4096*huffman_tree_num+1;
        auto huffman_bs_buffer=make_unique_for_overwrite<unsigned char[]>(huffman_bs_capacity);
        unsigned char *huffman_bs_begin=huffman_bs_buffer.get(),*huffman_bs=huffman_bs_begin;
        *huffman_bs++=selected_plan_id|0xc0;
        if constexpr(measure_phase_time)
            phase_start_time=omp_get_wtime();
        for(uint32_t tree_id=0;tree_id<huffman_tree_num;tree_id++)
        {
            int32_t *tree_qb=block_qb+tree_offsets[tree_id];
            uint64_t tree_qb_num=tree_offsets[tree_id+1]-tree_offsets[tree_id];
            mrcz::huffman_encoder_class<int32_t> encoder;
            encoder.preprocess_encode(tree_qb,tree_qb_num);
            encoder.save(huffman_bs);
            encoder.encode(tree_qb,tree_qb_num,huffman_bs);
        }
        if(!unpredictable.empty())
        {
            memcpy(huffman_bs,unpredictable.data(),unpredictable.size()*sizeof(T));
            huffman_bs+=unpredictable.size()*sizeof(T);
        }
        if constexpr(measure_phase_time)
            huffman_time[block_id]=omp_get_wtime()-phase_start_time;
        uint64_t huffman_bs_len=huffman_bs-huffman_bs_begin;
        if constexpr(measure_phase_time)
            phase_start_time=omp_get_wtime();
        uint64_t max_zstd_bs_len=ZSTD_compressBound(huffman_bs_len);
        auto zstd_bs_buffer=make_unique_for_overwrite<unsigned char[]>(max_zstd_bs_len);
        uint64_t zstd_bs_len=ZSTD_compress(zstd_bs_buffer.get(),max_zstd_bs_len,huffman_bs_begin,huffman_bs_len,3);
        block.zstd_bs.assign(zstd_bs_buffer.get(),zstd_bs_buffer.get()+zstd_bs_len);
        if constexpr(measure_phase_time)
            zstd_time[block_id]=omp_get_wtime()-phase_start_time;
    }

    // 执行单个field的完整FDWZ压缩；模板参数决定是否编译分段计时代码。
    template<bool measure_phase_time,typename T>
    void fdwz_compress_impl(mrcz::param_class<T> &params,string field,uint32_t huffman_tree_num,string fixed_plan)
    {
        auto &header=params.header;
        auto &blocks=params.partitions;
        auto &data=params.data;
        auto &x0=params.x[0];
        auto &x1=params.x[1];
        auto &x2=params.x[2];

        uint32_t n0=header.c_element_shape[0];
        uint32_t n1=header.c_element_shape[1];
        uint32_t n2=header.c_element_shape[2];
        uint32_t n12=n1*n2;
        uint32_t n012=n0*n1*n2;
        uint32_t corner_num=1;
        for (uint32_t dim=0;dim<3;dim++)
            if (header.c_element_shape[dim]>1)
                corner_num*=2;
        uint32_t residual_num=n012-corner_num;
        T abs_eb=params.abs_eb[field];
        T neg_2_eb=1/(2*abs_eb);
        long double field_max_abs=max(abs(static_cast<long double>(params.data_range_min[field])),abs(static_cast<long double>(params.data_range_max[field])));
        bool strict_error_check=field_max_abs*numeric_limits<T>::epsilon()*64>static_cast<long double>(abs_eb);
        // 与SZ3计时口径一致，总压缩时间包含内部工作区分配、plan生成、决策、预测量化和无损压缩。
        double start_time=omp_get_wtime();

        // 分别建立三个坐标方向的物理中点hierarchy。长度为1的方向只有一个始终可用的坐标点。
        vector<vector<uint32_t>> left_indices0,right_indices0,mid_indices0,compressed_indices0;
        vector<vector<uint32_t>> left_indices1,right_indices1,mid_indices1,compressed_indices1;
        vector<vector<uint32_t>> left_indices2,right_indices2,mid_indices2,compressed_indices2;
        left_indices0.push_back({0});       left_indices1.push_back({0});       left_indices2.push_back({0});
        right_indices0.push_back({n0-1});   right_indices1.push_back({n1-1});   right_indices2.push_back({n2-1});
        if (n0==1)
            compressed_indices0.push_back({0});
        else
        {
            compressed_indices0.push_back({0,n0-1});
            fdwz_generate_hierarchy(left_indices0,right_indices0,mid_indices0,compressed_indices0,x0);
        }
        if (n1==1)
            compressed_indices1.push_back({0});
        else
        {
            compressed_indices1.push_back({0,n1-1});
            fdwz_generate_hierarchy(left_indices1,right_indices1,mid_indices1,compressed_indices1,x1);
        }
        if (n2==1)
            compressed_indices2.push_back({0});
        else
        {
            compressed_indices2.push_back({0,n2-1});
            fdwz_generate_hierarchy(left_indices2,right_indices2,mid_indices2,compressed_indices2,x2);
        }
        uint32_t max_level=max(max(mid_indices0.size(),mid_indices1.size()),mid_indices2.size());
        vector<uint32_t> corner_positions=fdwz_generate_corner_positions(n0,n1,n2);

        // 预先构造6种维度顺序乘以linear/cubic两种插值方式，共12套可复用预测结构。
        vector<fdwz_prediction_plan<T>> sampled_plans;
        auto all_orders=fdwz_all_dimension_orders();
        for (uint32_t cubic_id=0;cubic_id<2;cubic_id++)
        {
            bool use_cubic=cubic_id==1;
            auto interpolation_stencils0=fdwz_generate_interpolation_stencils(left_indices0,right_indices0,mid_indices0,compressed_indices0,x0,use_cubic);
            auto interpolation_stencils1=fdwz_generate_interpolation_stencils(left_indices1,right_indices1,mid_indices1,compressed_indices1,x1,use_cubic);
            auto interpolation_stencils2=fdwz_generate_interpolation_stencils(left_indices2,right_indices2,mid_indices2,compressed_indices2,x2,use_cubic);
            for (auto &order:all_orders)
            {
                fdwz_prediction_plan<T> plan;
                plan.stencils=fdwz_generate_prediction_stencils(n0,n1,n2,n12,max_level,mid_indices0,compressed_indices0,interpolation_stencils0,
                    mid_indices1,compressed_indices1,interpolation_stencils1,mid_indices2,compressed_indices2,interpolation_stencils2,plan.level_offsets,order);
                if(!sampled_plans.empty() and plan.level_offsets!=sampled_plans[0].level_offsets)
                    throw runtime_error("FDWZ采样配置的预测层结构不一致");
                sampled_plans.push_back(move(plan));
            }
        }
        bool residual_qb_range_safe=fdwz_residual_qb_range_is_safe(params.data_range_min[field],params.data_range_max[field],abs_eb,sampled_plans);
        strict_error_check=strict_error_check or !residual_qb_range_safe;
        int32_t fixed_plan_id=-1;
        if (fixed_plan!="sample")
            for (uint32_t plan_id=0;plan_id<sampled_plans.size();plan_id++)
            {
                string plan_name;
                for (auto dimension:all_orders[plan_id%6])
                    plan_name.push_back("XYZ"[dimension]);
                plan_name+=plan_id<6?"-linear":"-cubic";
                if (plan_name==fixed_plan)
                    fixed_plan_id=plan_id;
            }
        if (fixed_plan_id<0 and fixed_plan!="sample")
            throw runtime_error("未知的FDWZ固定预测方案: "+fixed_plan);
        auto &prediction_level_offsets=sampled_plans[0].level_offsets;
        // 0表示正式版自动选择最多5棵tree；消融版可以强制只使用一棵tree。
        if (huffman_tree_num==0)
            huffman_tree_num=min<uint32_t>(5,prediction_level_offsets.size());
        if (huffman_tree_num<1 or huffman_tree_num>5)
            throw runtime_error("FDWZ Huffman tree数量必须为1、2、3、4或5");
        vector<uint32_t> residual_group_offsets;
        if (huffman_tree_num>=2)
        {
            if (prediction_level_offsets.size()<huffman_tree_num)
                throw runtime_error("FDWZ Huffman tree数量超过当前element的角点和预测层总数");
            residual_group_offsets.push_back(0);
            for (uint32_t level=1;level<huffman_tree_num-1;level++)
                residual_group_offsets.push_back(prediction_level_offsets[level]);
            residual_group_offsets.push_back(residual_num);
        }
        double structure_generation_time=0;
        if constexpr(measure_phase_time)
            structure_generation_time=omp_get_wtime()-start_time;
        vector<unique_ptr<int32_t[]>> partition_qb(blocks.size());
        for(uint64_t block_id=0;block_id<blocks.size();block_id++)
            partition_qb[block_id]=make_unique_for_overwrite<int32_t[]>(n012*blocks[block_id].element_num);
        vector<unsigned char> selected_plan_ids(blocks.size(),0);
        vector<uint64_t> sample_element_counts(blocks.size(),0);
        vector<uint64_t> unpredictable_counts(blocks.size(),0);
        array<uint32_t,5> tree_value_num{};
        if (huffman_tree_num==1)
            tree_value_num[0]=n012;
        else
        {
            tree_value_num[0]=corner_num;
            for (uint32_t group_id=0;group_id+1<residual_group_offsets.size();group_id++)
                tree_value_num[group_id+1]=residual_group_offsets[group_id+1]-residual_group_offsets[group_id];
        }
        vector<double> decision_time(measure_phase_time?blocks.size():0,0),prediction_quantization_time(measure_phase_time?blocks.size():0,0),huffman_time(measure_phase_time?blocks.size():0,0),zstd_time(measure_phase_time?blocks.size():0,0);
        double decision_start_time=0,decision_wall_time=0;
        if constexpr(measure_phase_time)
            decision_start_time=omp_get_wtime();
        // decision和正式压缩共用同一个OpenMP线程团队；两个omp for之间的barrier保证全部plan选择完成后才开始写正式bitstream。
        #pragma omp parallel num_threads(params.max_cpus)
        {
            fdwz_sampling_workspace<T> sampling_workspace;
            // 每个partition独立采样，因此可以并行为不同空间区域选择不同plan；同一线程复用全部候选评估工作区。
            #pragma omp for schedule(static)
            for (uint64_t block_id=0;block_id<blocks.size();block_id++)
            {
                double block_decision_start_time=0;
                if constexpr(measure_phase_time)
                    block_decision_start_time=omp_get_wtime();
                if (fixed_plan_id>=0)
                    selected_plan_ids[block_id]=fixed_plan_id;
                else
                {
                    auto &block=blocks[block_id];
                    uint64_t element_num=block.end_index-block.begin_index;
                    uint64_t sample_num;
                    if(params.fdwz_sample_interval==0)
                        sample_num=min<uint64_t>(element_num,min<uint64_t>(1024,max<uint64_t>(128,element_num/4)));
                    else
                        sample_num=min<uint64_t>(element_num,max<uint64_t>(params.fdwz_min_sample_num,1+(element_num-1)/params.fdwz_sample_interval));
                    sample_element_counts[block_id]=sample_num;
                    auto &original_sample=sampling_workspace.original_sample;
                    original_sample.resize(sample_num*n012);
                    for (uint64_t sample_id=0;sample_id<sample_num;sample_id++)
                    {
                        uint64_t relative_element_id;
                        relative_element_id=((2*sample_id+1)*element_num)/(2*sample_num);
                        memcpy(original_sample.data()+sample_id*n012,data.data()+(block.begin_index+relative_element_id)*n012,n012*sizeof(T));
                    }
                    bool sample_linear=params.fdwz_selector!="cubic-exhaustive";
                    uint32_t candidate_begin=sample_linear?0:6;
                    uint32_t candidate_end=sample_linear?6:12;
                    uint32_t best_plan_id=candidate_begin;
                    uint64_t best_size=fdwz_sample_compressed_size(original_sample,n012,corner_positions,sampled_plans[best_plan_id],residual_group_offsets,huffman_tree_num,abs_eb,best_plan_id,sampling_workspace);
                    for(uint32_t plan_id=candidate_begin+1;plan_id<candidate_end;plan_id++)
                    {
                        uint64_t candidate_size=fdwz_sample_compressed_size(original_sample,n012,corner_positions,sampled_plans[plan_id],residual_group_offsets,huffman_tree_num,abs_eb,plan_id,sampling_workspace);
                        if(candidate_size<best_size)
                            best_size=candidate_size,best_plan_id=plan_id;
                    }
                    selected_plan_ids[block_id]=params.fdwz_selector=="linear-proxy"?best_plan_id+6:best_plan_id;
                }
                if constexpr(measure_phase_time)
                    decision_time[block_id]=omp_get_wtime()-block_decision_start_time;
            }
            #pragma omp master
            {
                if constexpr(measure_phase_time)
                    decision_wall_time=omp_get_wtime()-decision_start_time;
            }
            #pragma omp barrier
            // partition覆盖互不重叠的element和QB范围，因此各线程可以安全地原地覆盖data并写入各自的bitstream。
            #pragma omp for schedule(static)
            for (uint64_t block_id=0;block_id<blocks.size();block_id++)
            {
                auto &prediction_stencils=sampled_plans[selected_plan_ids[block_id]].stencils;
                if(strict_error_check)
                    fdwz_compress_block<measure_phase_time,true,true>(params,block_id,n012,corner_num,residual_num,corner_positions,prediction_stencils,residual_group_offsets,tree_value_num,huffman_tree_num,abs_eb,neg_2_eb,selected_plan_ids[block_id],partition_qb[block_id].get(),unpredictable_counts,prediction_quantization_time,huffman_time,zstd_time);
                else if(residual_qb_range_safe)
                    fdwz_compress_block<measure_phase_time,false,false>(params,block_id,n012,corner_num,residual_num,corner_positions,prediction_stencils,residual_group_offsets,tree_value_num,huffman_tree_num,abs_eb,neg_2_eb,selected_plan_ids[block_id],partition_qb[block_id].get(),unpredictable_counts,prediction_quantization_time,huffman_time,zstd_time);
                else
                    fdwz_compress_block<measure_phase_time,true,false>(params,block_id,n012,corner_num,residual_num,corner_positions,prediction_stencils,residual_group_offsets,tree_value_num,huffman_tree_num,abs_eb,neg_2_eb,selected_plan_ids[block_id],partition_qb[block_id].get(),unpredictable_counts,prediction_quantization_time,huffman_time,zstd_time);
            }
        }
        double end_time=omp_get_wtime();
        vector<uint64_t> configuration_counts(sampled_plans.size(),0);
        for (auto plan_id:selected_plan_ids)
            configuration_counts[plan_id]++;
        printf("Sampled configurations:");
        for (uint32_t plan_id=0;plan_id<configuration_counts.size();plan_id++)
            if (configuration_counts[plan_id]>0)
                printf(" %c%c%c-%s=%lu","XYZ"[all_orders[plan_id%6][0]],"XYZ"[all_orders[plan_id%6][1]],"XYZ"[all_orders[plan_id%6][2]],plan_id<6?"linear":"cubic",configuration_counts[plan_id]);
        printf("\n");
        printf("FDWZ QB width= 32 bits, unpredictable value num= %lu, error check= %s, residual range check= %s\n",accumulate(unpredictable_counts.begin(),unpredictable_counts.end(),uint64_t(0)),strict_error_check?"strict":"fast",residual_qb_range_safe?"field":"per-value");
        if(fixed_plan_id<0)
        {
            printf("FDWZ selector= %s\n",params.fdwz_selector.c_str());
            uint64_t total_sample_element_num=accumulate(sample_element_counts.begin(),sample_element_counts.end(),uint64_t(0));
            printf("FDWZ sampling interval= %lu, minimum samples per partition= %lu, sampled element num= %lu\n",params.fdwz_sample_interval,params.fdwz_min_sample_num,total_sample_element_num);
        }
        if constexpr(measure_phase_time)
        {
            // 分段时间是所有partition的CPU时间之和；Compression time则是包含并行加速后的wall time。
            double total_decision_time=0,total_prediction_quantization_time=0,total_huffman_time=0,total_zstd_time=0;
            for(uint64_t block_id=0;block_id<blocks.size();block_id++)
            {
                total_decision_time+=decision_time[block_id];
                total_prediction_quantization_time+=prediction_quantization_time[block_id];
                total_huffman_time+=huffman_time[block_id];
                total_zstd_time+=zstd_time[block_id];
            }
            printf("FDWZ prediction structure generation time= %f seconds\n",structure_generation_time);
            printf("FDWZ decision wall time= %f seconds, accumulated partition CPU time= %f seconds\n",decision_wall_time,total_decision_time);
            double measured_stage_time=total_decision_time+total_prediction_quantization_time+total_huffman_time+total_zstd_time;
            printf("FDWZ compression accumulated partition CPU times: prediction+quantization+QB placement= %f, Huffman= %f, Zstd= %f seconds\n",total_prediction_quantization_time,total_huffman_time,total_zstd_time);
            if(measured_stage_time>0)
                printf("FDWZ accumulated CPU-time shares: decision= %.2f%%, prediction+quantization+QB placement= %.2f%%, Huffman= %.2f%%, Zstd= %.2f%%\n",100*total_decision_time/measured_stage_time,100*total_prediction_quantization_time/measured_stage_time,100*total_huffman_time/measured_stage_time,100*total_zstd_time/measured_stage_time);
        }
        printf("Compression time= %f seconds\n",end_time-start_time);
        printf("\033[91mField= %s, compression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),params.element_num_to_process*n012*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
    }
}

namespace mrcz
{
    // 根据正式版或消融版配置选择Huffman tree数量，并进入FDWZ压缩实现。
    template<typename T>
    template<bool measure_phase_time,bool ablation>
    void param_class<T>::compress_fdwz(string field)
    {
        // 正式FDWZ忽略消融开关并使用自动tree数量；FDWZ-ablation才读取if_max_huffman_tree。
        uint32_t huffman_tree_num=0;
        if constexpr(ablation)
            huffman_tree_num=if_max_huffman_tree?0:1;
        string selected_plan="sample";
        if constexpr(ablation)
            selected_plan=fdwz_plan;
        detail::fdwz_compress_impl<measure_phase_time>(*this,field,huffman_tree_num,selected_plan);
    }
}
