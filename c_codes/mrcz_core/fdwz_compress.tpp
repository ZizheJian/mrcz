#pragma once


#if defined(__AVX2__)
#include <immintrin.h>
#endif
#include <set>
#include "mrcz.hpp"
#include "huffman.hpp"
#include "quantize.tpp"
#include "sz3_index_order_prediction.tpp"
#include "sz3_normal_boundary_prediction.tpp"
#include "sz3_surface_first_prediction.tpp"
#include <zstd.h>

template<typename T>
void sz3_generate_weight(vector<vector<uint32_t>> &left_indices,vector<vector<uint32_t>> &right_indices,vector<vector<uint32_t>> &mid_indices,
    vector<vector<uint32_t>> &compressed_indices,vector<vector<T>> &left_weight,vector<vector<T>> &right_weight,vector<T> &x,bool use_index_midpoint=false)
{
    while (true)
    {
        mid_indices.push_back({});
        left_weight.push_back({});
        right_weight.push_back({});
        vector<uint32_t> new_left_indices,new_right_indices;
        for (uint32_t i=0;i<left_indices.back().size();i++)
        {
            uint32_t left_index=left_indices.back()[i],right_index=right_indices.back()[i],mid_index=left_index+1;
            if (use_index_midpoint)
                mid_index=(left_index+right_index)/2;
            else
            {
                T ideal_mid_pos=(x[left_index]+x[right_index])/2;
                for (uint32_t j=left_index+2;j<right_index;j++)
                    if (fabs(x[j]-ideal_mid_pos)<fabs(x[mid_index]-ideal_mid_pos))
                        mid_index=j;
            }
            mid_indices.back().push_back(mid_index);
            left_weight.back().push_back((x[right_index]-x[mid_index])/(x[right_index]-x[left_index]));
            right_weight.back().push_back((x[mid_index]-x[left_index])/(x[right_index]-x[left_index]));
            if (mid_index-left_index>1)
                new_left_indices.push_back(left_index),new_right_indices.push_back(mid_index);
            if (right_index-mid_index>1)
                new_left_indices.push_back(mid_index),new_right_indices.push_back(right_index);
        }
        set<uint32_t> temp_set;
        temp_set.insert(compressed_indices.back().begin(),compressed_indices.back().end());
        temp_set.insert(mid_indices.back().begin(),mid_indices.back().end());
        compressed_indices.push_back(vector<uint32_t>(temp_set.begin(),temp_set.end()));
        if (new_left_indices.size()==0)
            break;
        else
            left_indices.push_back(new_left_indices),right_indices.push_back(new_right_indices);
    }
}

template<typename T>
struct sz3_cubic_stencil
{
    uint32_t point_num=0;
    uint32_t indices[4]={0,0,0,0};
    T weights[4]={0,0,0,0};
};

enum class sz3_avx2_rounding:uint8_t
{
    disabled,
    ties_to_even,
    ties_away_from_zero
};

template<typename T>
vector<vector<sz3_cubic_stencil<T>>> sz3_generate_cubic_stencils(vector<vector<uint32_t>> &left_indices,vector<vector<uint32_t>> &right_indices,
    vector<vector<uint32_t>> &mid_indices,vector<vector<uint32_t>> &compressed_indices,vector<T> &x)
{
    vector<vector<sz3_cubic_stencil<T>>> stencils(left_indices.size());
    for (uint32_t level=0;level<left_indices.size();level++)
    {
        stencils[level].resize(left_indices[level].size());
        if (compressed_indices[level].size()<3)
            continue;
        for (uint32_t interval=0;interval<left_indices[level].size();interval++)
        {
            auto &stencil=stencils[level][interval];
            uint32_t left_index=left_indices[level][interval],right_index=right_indices[level][interval],mid_index=mid_indices[level][interval];
            uint32_t outer_left=left_index,outer_right=right_index;
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

inline vector<uint32_t> sz3_dimension_order_from_string(string order)
{
    if (order.size()!=3)
        throw runtime_error("SZ3维度顺序必须由X、Y、Z各出现一次组成");
    vector<uint32_t> dimension_order;
    vector<bool> seen(3,false);
    for (auto dimension_name:order)
    {
        uint32_t dimension=dimension_name=='X'?0:dimension_name=='Y'?1:dimension_name=='Z'?2:3;
        if (dimension>=3 or seen[dimension])
            throw runtime_error("SZ3维度顺序必须由X、Y、Z各出现一次组成");
        seen[dimension]=true;
        dimension_order.push_back(dimension);
    }
    return dimension_order;
}

inline vector<uint32_t> sz3_dimension_order_from_method(string method)
{
    string prefix="FDWZ-order-";
    if (method.rfind(prefix,0)!=0)
        return {0,1,2};
    return sz3_dimension_order_from_string(method.substr(prefix.size()));
}

struct sz3_sampled_configuration
{
    vector<uint32_t> dimension_order;
    bool use_cubic;
    uint32_t huffman_tree_num;
};

inline sz3_sampled_configuration sz3_configuration_from_suffix(string suffix)
{
    size_t first_separator=suffix.find('-');
    size_t second_separator=first_separator==string::npos?string::npos:suffix.find('-',first_separator+1);
    if (first_separator==string::npos or second_separator==string::npos or suffix.find('-',second_separator+1)!=string::npos)
        throw runtime_error("FDWZ采样配置方法名应为FDWZ-sampled-<XYZ顺序>-<linear|cubic>-<1..5>");
    string order=suffix.substr(0,first_separator);
    string predictor=suffix.substr(first_separator+1,second_separator-first_separator-1);
    string tree_num_string=suffix.substr(second_separator+1);
    if (predictor!="linear" and predictor!="cubic")
        throw runtime_error("SZ3采样配置的预测函数必须为linear或cubic");
    uint32_t tree_num=0;
    auto parse_result=from_chars(tree_num_string.data(),tree_num_string.data()+tree_num_string.size(),tree_num);
    if (parse_result.ec!=errc() or parse_result.ptr!=tree_num_string.data()+tree_num_string.size() or tree_num<1 or tree_num>5)
        throw runtime_error("SZ3采样配置的Huffman tree数量必须为1到5");
    return {sz3_dimension_order_from_string(order),predictor=="cubic",tree_num};
}

inline sz3_sampled_configuration sz3_sampled_configuration_from_method(string method)
{
    string prefix="FDWZ-sampled-";
    if (method.rfind(prefix,0)!=0)
        throw runtime_error("FDWZ采样配置方法名必须以FDWZ-sampled-开头");
    return sz3_configuration_from_suffix(method.substr(prefix.size()));
}

inline sz3_sampled_configuration sz3_position_bias_configuration_from_method(string method)
{
    string prefix="FDWZ-position-bias-";
    if (method.rfind(prefix,0)!=0)
        throw runtime_error("FDWZ位置bias方法名必须以FDWZ-position-bias-开头");
    return sz3_configuration_from_suffix(method.substr(prefix.size()));
}

inline vector<uint32_t> sz3_generate_prediction_positions(uint32_t n0,uint32_t n1,uint32_t n2,uint32_t n12,uint32_t max_level,
    vector<vector<uint32_t>> &left_indices0,vector<vector<uint32_t>> &right_indices0,vector<vector<uint32_t>> &mid_indices0,vector<vector<uint32_t>> &compressed_indices0,
    vector<vector<uint32_t>> &left_indices1,vector<vector<uint32_t>> &right_indices1,vector<vector<uint32_t>> &mid_indices1,vector<vector<uint32_t>> &compressed_indices1,
    vector<vector<uint32_t>> &left_indices2,vector<vector<uint32_t>> &right_indices2,vector<vector<uint32_t>> &mid_indices2,vector<vector<uint32_t>> &compressed_indices2,
    vector<uint32_t> &level_offsets,vector<uint32_t> dimension_order,vector<bool> *reference_mask=nullptr)
{
    vector<uint32_t> positions;
    level_offsets={0};
    if (reference_mask!=nullptr)
        reference_mask->assign(n0*n1*n2,false);
    uint32_t level0=0,level1=0,level2=0;
    for (uint32_t level=0;level<max_level;level++)
    {
        for (auto dimension:dimension_order)
        {
            if (dimension==0 and n0>1 and level0<left_indices0.size())
            {
                for (uint32_t i0=0;i0<left_indices0[level0].size();i0++)
                    for (auto i1:compressed_indices1[level1])
                        for (auto i2:compressed_indices2[level2])
                        {
                            positions.push_back(i1*n2+i2+mid_indices0[level0][i0]*n12);
                            if (reference_mask!=nullptr)
                            {
                                (*reference_mask)[i1*n2+i2+left_indices0[level0][i0]*n12]=true;
                                (*reference_mask)[i1*n2+i2+right_indices0[level0][i0]*n12]=true;
                            }
                        }
                level0++;
            }
            else if (dimension==1 and n1>1 and level1<left_indices1.size())
            {
                for (auto i0:compressed_indices0[level0])
                    for (uint32_t i1=0;i1<left_indices1[level1].size();i1++)
                        for (auto i2:compressed_indices2[level2])
                        {
                            positions.push_back(i0*n12+i2+mid_indices1[level1][i1]*n2);
                            if (reference_mask!=nullptr)
                            {
                                (*reference_mask)[i0*n12+i2+left_indices1[level1][i1]*n2]=true;
                                (*reference_mask)[i0*n12+i2+right_indices1[level1][i1]*n2]=true;
                            }
                        }
                level1++;
            }
            else if (dimension==2 and n2>1 and level2<left_indices2.size())
            {
                for (auto i0:compressed_indices0[level0])
                    for (auto i1:compressed_indices1[level1])
                        for (uint32_t i2=0;i2<left_indices2[level2].size();i2++)
                        {
                            positions.push_back(i0*n12+i1*n2+mid_indices2[level2][i2]);
                            if (reference_mask!=nullptr)
                            {
                                (*reference_mask)[i0*n12+i1*n2+left_indices2[level2][i2]]=true;
                                (*reference_mask)[i0*n12+i1*n2+right_indices2[level2][i2]]=true;
                            }
                        }
                level2++;
            }
        }
        level_offsets.push_back(positions.size());
    }
    uint32_t corner_num=1;
    for (uint32_t dim_size:{n0,n1,n2})
        if (dim_size>1)
            corner_num*=2;
    uint32_t n012=n0*n1*n2;
    if (positions.size()!=n012-corner_num)
        throw runtime_error("SZ3 prediction position数量错误");
    vector<bool> seen(n012,false);
    for (auto position:positions)
    {
        if (position>=n012 or seen[position])
            throw runtime_error("SZ3 prediction position重复或越界");
        seen[position]=true;
    }
    return positions;
}

inline vector<uint32_t> sz3_generate_corner_positions(uint32_t n0,uint32_t n1,uint32_t n2)
{
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

template<typename T>
struct sz3_prediction_stencil
{
    uint32_t target=0;
    uint32_t point_num=0;
    uint32_t indices[4]={0,0,0,0};
    T weights[4]={0,0,0,0};
};

template<typename T>
inline T sz3_predict_from_stencil(T *element,sz3_prediction_stencil<T> &stencil)
{
    T prediction=0;
    for (uint32_t i=0;i<stencil.point_num;i++)
        prediction+=stencil.weights[i]*element[stencil.indices[i]];
    return prediction;
}

template<sz3_avx2_rounding rounding_mode,typename T>
inline void sz3_quantize_scalar(T pred,T &truth,int64_t &qb,T abs_eb,T neg_2_eb)
{
    T scaled=(pred-truth)*neg_2_eb;
    if constexpr (rounding_mode==sz3_avx2_rounding::ties_to_even)
        qb=static_cast<int64_t>(nearbyint(scaled));
    else
        qb=static_cast<int64_t>(round(scaled));
    truth=pred-qb*2*abs_eb;
}

#if defined(__AVX2__)
template<sz3_avx2_rounding rounding_mode>
inline bool sz3_quantize_stencil_avx2(float *first_element,uint32_t n012,sz3_prediction_stencil<float> &stencil,float abs_eb,float neg_2_eb,int64_t *qb,float *reconstructed)
{
    __m256i offsets=_mm256_setr_epi32(0,n012,2*n012,3*n012,4*n012,5*n012,6*n012,7*n012);
    __m256 prediction=_mm256_setzero_ps();
    for (uint32_t reference=0;reference<stencil.point_num;reference++)
    {
        __m256 values=_mm256_i32gather_ps(first_element+stencil.indices[reference],offsets,4);
        prediction=_mm256_add_ps(prediction,_mm256_mul_ps(_mm256_set1_ps(stencil.weights[reference]),values));
    }
    __m256 truth=_mm256_i32gather_ps(first_element+stencil.target,offsets,4);
    __m256 scaled=_mm256_mul_ps(_mm256_sub_ps(prediction,truth),_mm256_set1_ps(neg_2_eb));
    __m256 sign_mask=_mm256_set1_ps(-0.0f);
    __m256 magnitude=_mm256_andnot_ps(sign_mask,scaled);
    __m256 in_range=_mm256_cmp_ps(magnitude,_mm256_set1_ps(2147483520.0f),_CMP_LE_OQ);
    if (_mm256_movemask_ps(in_range)!=0xff)
        return false;
    __m256i qb32;
    if constexpr (rounding_mode==sz3_avx2_rounding::ties_to_even)
        qb32=_mm256_cvttps_epi32(_mm256_round_ps(scaled,_MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC));
    else
    {
        __m256 half_with_sign=_mm256_or_ps(_mm256_and_ps(scaled,sign_mask),_mm256_set1_ps(0.5f));
        qb32=_mm256_cvttps_epi32(_mm256_add_ps(scaled,half_with_sign));
    }
    __m256 reconstructed_vector=_mm256_sub_ps(prediction,_mm256_mul_ps(_mm256_cvtepi32_ps(qb32),_mm256_set1_ps(2*abs_eb)));
    alignas(32) int32_t qb32_array[8];
    _mm256_store_si256(reinterpret_cast<__m256i*>(qb32_array),qb32);
    _mm256_store_ps(reconstructed,reconstructed_vector);
    for (uint32_t lane=0;lane<8;lane++)
        qb[lane]=qb32_array[lane];
    return true;
}

template<sz3_avx2_rounding rounding_mode>
inline bool sz3_quantize_stencil_avx2(double *first_element,uint32_t n012,sz3_prediction_stencil<double> &stencil,double abs_eb,double neg_2_eb,int64_t *qb,double *reconstructed)
{
    __m128i offsets=_mm_setr_epi32(0,n012,2*n012,3*n012);
    __m256d prediction=_mm256_setzero_pd();
    for (uint32_t reference=0;reference<stencil.point_num;reference++)
    {
        __m256d values=_mm256_i32gather_pd(first_element+stencil.indices[reference],offsets,8);
        prediction=_mm256_add_pd(prediction,_mm256_mul_pd(_mm256_set1_pd(stencil.weights[reference]),values));
    }
    __m256d truth=_mm256_i32gather_pd(first_element+stencil.target,offsets,8);
    __m256d scaled=_mm256_mul_pd(_mm256_sub_pd(prediction,truth),_mm256_set1_pd(neg_2_eb));
    __m256d sign_mask=_mm256_set1_pd(-0.0);
    __m256d magnitude=_mm256_andnot_pd(sign_mask,scaled);
    __m256d in_range=_mm256_cmp_pd(magnitude,_mm256_set1_pd(2147483647.0),_CMP_LE_OQ);
    if (_mm256_movemask_pd(in_range)!=0xf)
        return false;
    __m128i qb32;
    if constexpr (rounding_mode==sz3_avx2_rounding::ties_to_even)
        qb32=_mm256_cvttpd_epi32(_mm256_round_pd(scaled,_MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC));
    else
    {
        __m256d half_with_sign=_mm256_or_pd(_mm256_and_pd(scaled,sign_mask),_mm256_set1_pd(0.5));
        qb32=_mm256_cvttpd_epi32(_mm256_add_pd(scaled,half_with_sign));
    }
    __m256d reconstructed_vector=_mm256_sub_pd(prediction,_mm256_mul_pd(_mm256_cvtepi32_pd(qb32),_mm256_set1_pd(2*abs_eb)));
    alignas(16) int32_t qb32_array[4];
    _mm_store_si128(reinterpret_cast<__m128i*>(qb32_array),qb32);
    _mm256_store_pd(reconstructed,reconstructed_vector);
    for (uint32_t lane=0;lane<4;lane++)
        qb[lane]=qb32_array[lane];
    return true;
}
#endif

template<typename T>
vector<sz3_prediction_stencil<T>> sz3_generate_prediction_stencils(uint32_t n0,uint32_t n1,uint32_t n2,uint32_t n12,uint32_t max_level,
    vector<vector<uint32_t>> &left_indices0,vector<vector<uint32_t>> &right_indices0,vector<vector<uint32_t>> &mid_indices0,vector<vector<uint32_t>> &compressed_indices0,vector<vector<T>> &left_weight0,vector<vector<T>> &right_weight0,vector<vector<sz3_cubic_stencil<T>>> &cubic_stencils0,
    vector<vector<uint32_t>> &left_indices1,vector<vector<uint32_t>> &right_indices1,vector<vector<uint32_t>> &mid_indices1,vector<vector<uint32_t>> &compressed_indices1,vector<vector<T>> &left_weight1,vector<vector<T>> &right_weight1,vector<vector<sz3_cubic_stencil<T>>> &cubic_stencils1,
    vector<vector<uint32_t>> &left_indices2,vector<vector<uint32_t>> &right_indices2,vector<vector<uint32_t>> &mid_indices2,vector<vector<uint32_t>> &compressed_indices2,vector<vector<T>> &left_weight2,vector<vector<T>> &right_weight2,vector<vector<sz3_cubic_stencil<T>>> &cubic_stencils2,
    vector<uint32_t> dimension_order,bool use_cubic)
{
    vector<sz3_prediction_stencil<T>> stencils;
    uint32_t level0=0,level1=0,level2=0;
    for (uint32_t level=0;level<max_level;level++)
        for (auto dimension:dimension_order)
        {
            if (dimension==0 and n0>1 and level0<left_indices0.size())
            {
                for (uint32_t interval=0;interval<left_indices0[level0].size();interval++)
                    for (auto i1:compressed_indices1[level1])
                        for (auto i2:compressed_indices2[level2])
                        {
                            uint32_t pos_base=i1*n2+i2;
                            sz3_prediction_stencil<T> stencil;
                            stencil.target=pos_base+mid_indices0[level0][interval]*n12;
                            if (use_cubic and cubic_stencils0[level0][interval].point_num>0)
                            {
                                stencil.point_num=cubic_stencils0[level0][interval].point_num;
                                for (uint32_t i=0;i<stencil.point_num;i++)
                                    stencil.indices[i]=pos_base+cubic_stencils0[level0][interval].indices[i]*n12,stencil.weights[i]=cubic_stencils0[level0][interval].weights[i];
                            }
                            else
                            {
                                stencil.point_num=2;
                                stencil.indices[0]=pos_base+left_indices0[level0][interval]*n12;
                                stencil.indices[1]=pos_base+right_indices0[level0][interval]*n12;
                                stencil.weights[0]=left_weight0[level0][interval];
                                stencil.weights[1]=right_weight0[level0][interval];
                            }
                            stencils.push_back(stencil);
                        }
                level0++;
            }
            else if (dimension==1 and n1>1 and level1<left_indices1.size())
            {
                for (auto i0:compressed_indices0[level0])
                    for (uint32_t interval=0;interval<left_indices1[level1].size();interval++)
                        for (auto i2:compressed_indices2[level2])
                        {
                            uint32_t pos_base=i0*n12+i2;
                            sz3_prediction_stencil<T> stencil;
                            stencil.target=pos_base+mid_indices1[level1][interval]*n2;
                            if (use_cubic and cubic_stencils1[level1][interval].point_num>0)
                            {
                                stencil.point_num=cubic_stencils1[level1][interval].point_num;
                                for (uint32_t i=0;i<stencil.point_num;i++)
                                    stencil.indices[i]=pos_base+cubic_stencils1[level1][interval].indices[i]*n2,stencil.weights[i]=cubic_stencils1[level1][interval].weights[i];
                            }
                            else
                            {
                                stencil.point_num=2;
                                stencil.indices[0]=pos_base+left_indices1[level1][interval]*n2;
                                stencil.indices[1]=pos_base+right_indices1[level1][interval]*n2;
                                stencil.weights[0]=left_weight1[level1][interval];
                                stencil.weights[1]=right_weight1[level1][interval];
                            }
                            stencils.push_back(stencil);
                        }
                level1++;
            }
            else if (dimension==2 and n2>1 and level2<left_indices2.size())
            {
                for (auto i0:compressed_indices0[level0])
                    for (auto i1:compressed_indices1[level1])
                        for (uint32_t interval=0;interval<left_indices2[level2].size();interval++)
                        {
                            uint32_t pos_base=i0*n12+i1*n2;
                            sz3_prediction_stencil<T> stencil;
                            stencil.target=pos_base+mid_indices2[level2][interval];
                            if (use_cubic and cubic_stencils2[level2][interval].point_num>0)
                            {
                                stencil.point_num=cubic_stencils2[level2][interval].point_num;
                                for (uint32_t i=0;i<stencil.point_num;i++)
                                    stencil.indices[i]=pos_base+cubic_stencils2[level2][interval].indices[i],stencil.weights[i]=cubic_stencils2[level2][interval].weights[i];
                            }
                            else
                            {
                                stencil.point_num=2;
                                stencil.indices[0]=pos_base+left_indices2[level2][interval];
                                stencil.indices[1]=pos_base+right_indices2[level2][interval];
                                stencil.weights[0]=left_weight2[level2][interval];
                                stencil.weights[1]=right_weight2[level2][interval];
                            }
                            stencils.push_back(stencil);
                        }
                level2++;
            }
        }
    return stencils;
}

template<typename T>
struct sz3_prediction_plan
{
    vector<uint32_t> dimension_order;
    bool use_cubic=false;
    vector<uint32_t> positions;
    vector<uint32_t> level_offsets;
    vector<sz3_prediction_stencil<T>> stencils;
};

inline vector<vector<uint32_t>> sz3_all_dimension_orders()
{
    return {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
}

template<typename T>
uint64_t sz3_sample_compressed_size(vector<T> &original_sample,uint32_t n012,vector<uint32_t> &corner_positions,sz3_prediction_plan<T> &plan,
    vector<uint32_t> &residual_group_offsets,uint32_t huffman_tree_num,T abs_eb,unsigned char configuration)
{
    uint64_t sample_num=original_sample.size()/n012;
    uint32_t seed_num=corner_positions.size();
    uint32_t residual_num=n012-seed_num;
    T neg_2_eb=1/(2*abs_eb);
    vector<T> sample=original_sample;
    vector<int64_t> corner_qb(sample_num*seed_num),residual_qb(sample_num*residual_num);
    for (uint64_t sample_id=0;sample_id<sample_num;sample_id++)
    {
        T *element=sample.data()+sample_id*n012;
        for (uint32_t corner_id=0;corner_id<seed_num;corner_id++)
        {
            uint32_t target=corner_positions[corner_id];
            int64_t qb=static_cast<int64_t>(round(element[target]*neg_2_eb));
            corner_qb[sample_id*seed_num+corner_id]=qb;
            element[target]=qb*2*abs_eb;
        }
        for (uint32_t residual_id=0;residual_id<residual_num;residual_id++)
        {
            auto &stencil=plan.stencils[residual_id];
            T prediction=sz3_predict_from_stencil(element,stencil);
            int64_t qb=static_cast<int64_t>(round((prediction-element[stencil.target])*neg_2_eb));
            residual_qb[sample_id*residual_num+residual_id]=qb;
            element[stencil.target]=prediction-qb*2*abs_eb;
        }
    }
    vector<vector<int64_t>> tree_qb(huffman_tree_num);
    if (huffman_tree_num==1)
    {
        tree_qb[0].resize(sample_num*n012);
        for (uint64_t sample_id=0;sample_id<sample_num;sample_id++)
        {
            copy(corner_qb.begin()+sample_id*seed_num,corner_qb.begin()+(sample_id+1)*seed_num,tree_qb[0].begin()+sample_id*n012);
            copy(residual_qb.begin()+sample_id*residual_num,residual_qb.begin()+(sample_id+1)*residual_num,tree_qb[0].begin()+sample_id*n012+seed_num);
        }
    }
    else
    {
        tree_qb[0]=move(corner_qb);
        if (huffman_tree_num==2)
            tree_qb[1]=move(residual_qb);
        else
            for (uint32_t group_id=0;group_id+1<residual_group_offsets.size();group_id++)
            {
                uint32_t group_begin=residual_group_offsets[group_id],group_size=residual_group_offsets[group_id+1]-group_begin;
                tree_qb[group_id+1].resize(sample_num*group_size);
                for (uint64_t sample_id=0;sample_id<sample_num;sample_id++)
                    copy(residual_qb.begin()+sample_id*residual_num+group_begin,residual_qb.begin()+sample_id*residual_num+group_begin+group_size,tree_qb[group_id+1].begin()+sample_id*group_size);
            }
    }
    vector<unsigned char> huffman_bs_vec(sample_num*n012*20+4096*huffman_tree_num+1);
    unsigned char *huffman_bs=huffman_bs_vec.data();
    *huffman_bs++=configuration;
    for (auto &values:tree_qb)
    {
        mrcz::huffman_encoder_class<int64_t> encoder;
        encoder.preprocess_encode(values.data(),values.size());
        encoder.save(huffman_bs);
        encoder.encode(values.data(),values.size(),huffman_bs);
    }
    uint64_t huffman_bs_len=huffman_bs-huffman_bs_vec.data();
    vector<unsigned char> zstd_bs(ZSTD_compressBound(huffman_bs_len));
    uint64_t zstd_bs_len=ZSTD_compress(zstd_bs.data(),zstd_bs.size(),huffman_bs_vec.data(),huffman_bs_len,3);
    if (ZSTD_isError(zstd_bs_len))
        throw runtime_error("ZSTD_compress error while sampling SZ3 configuration");
    return zstd_bs_len;
}

inline double sz3_reference_eb_scale_from_method(string method)
{
    string prefix="FDWZ-reference-EB-";
    if (method.rfind(prefix,0)!=0)
        throw runtime_error("SZ3 reference EB方法名错误");
    string suffix=method.substr(prefix.size());
    size_t parsed_num=0;
    double scale=stod(suffix,&parsed_num);
    if (parsed_num!=suffix.size() or scale<=0 or scale>1)
        throw runtime_error("SZ3 reference EB比例必须大于0且不超过1");
    return scale;
}

template<typename T,sz3_avx2_rounding avx2_rounding=sz3_avx2_rounding::disabled>
void sz3_3d_gll_compress_impl(mrcz::param_class<T> &params,string field,uint32_t huffman_tree_num,vector<uint32_t> dimension_order={0,1,2},bool use_cubic=false,bool use_surface_first=false,bool use_normal_boundary=false,bool use_index_order=false,bool use_index_midpoint=false,T reference_eb_scale=1,bool use_position_bias=false,bool sample_configuration=false)
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
    uint32_t seed_num=use_index_order?1:corner_num;
    uint32_t residual_num=n012-seed_num;
    T abs_eb=params.abs_eb[field];
    T neg_2_eb=1/(2*abs_eb);
    vector<vector<uint32_t>> left_indices0,right_indices0,mid_indices0,compressed_indices0;
    vector<vector<uint32_t>> left_indices1,right_indices1,mid_indices1,compressed_indices1;
    vector<vector<uint32_t>> left_indices2,right_indices2,mid_indices2,compressed_indices2;
    vector<vector<T>> left_weight0,right_weight0;
    vector<vector<T>> left_weight1,right_weight1;
    vector<vector<T>> left_weight2,right_weight2;
    left_indices0.push_back({0});       left_indices1.push_back({0});       left_indices2.push_back({0});
    right_indices0.push_back({n0-1});   right_indices1.push_back({n1-1});   right_indices2.push_back({n2-1});
    if (n0==1)
        compressed_indices0.push_back({0});
    else
    {
        compressed_indices0.push_back({0,n0-1});
        sz3_generate_weight(left_indices0,right_indices0,mid_indices0,compressed_indices0,left_weight0,right_weight0,x0,use_index_midpoint);
    }
    if (n1==1)
        compressed_indices1.push_back({0});
    else
    {
        compressed_indices1.push_back({0,n1-1});
        sz3_generate_weight(left_indices1,right_indices1,mid_indices1,compressed_indices1,left_weight1,right_weight1,x1,use_index_midpoint);
    }
    if (n2==1)
        compressed_indices2.push_back({0});
    else
    {
        compressed_indices2.push_back({0,n2-1});
        sz3_generate_weight(left_indices2,right_indices2,mid_indices2,compressed_indices2,left_weight2,right_weight2,x2,use_index_midpoint);
    }
    uint32_t max_level=max(max(left_indices0.size(),left_indices1.size()),left_indices2.size());
    vector<vector<sz3_cubic_stencil<T>>> cubic_stencils0,cubic_stencils1,cubic_stencils2;
    if (use_cubic or sample_configuration)
    {
        cubic_stencils0=sz3_generate_cubic_stencils(left_indices0,right_indices0,mid_indices0,compressed_indices0,x0);
        cubic_stencils1=sz3_generate_cubic_stencils(left_indices1,right_indices1,mid_indices1,compressed_indices1,x1);
        cubic_stencils2=sz3_generate_cubic_stencils(left_indices2,right_indices2,mid_indices2,compressed_indices2,x2);
    }
    vector<uint32_t> prediction_level_offsets,prediction_positions;
    vector<bool> reference_mask;
    vector<sz3_index_order_stencil<T>> index_order_stencils;
    vector<sz3_surface_first_stencil<T>> surface_first_stencils;
    vector<sz3_normal_boundary_stencil<T>> normal_boundary_stencils;
    if (use_index_order)
    {
        index_order_stencils=sz3_generate_index_order_stencils(header.c_element_shape,params.x);
        for (auto &stencil:index_order_stencils)
            prediction_positions.push_back(stencil.target);
    }
    else if (use_normal_boundary)
    {
        normal_boundary_stencils=sz3_generate_normal_boundary_stencils(header.c_element_shape,params.x);
        for (auto &stencil:normal_boundary_stencils)
            prediction_positions.push_back(stencil.target);
    }
    else if (use_surface_first)
    {
        surface_first_stencils=sz3_generate_surface_first_stencils(header.c_element_shape,params.x);
        for (auto &stencil:surface_first_stencils)
            prediction_positions.push_back(stencil.target);
    }
    else
        prediction_positions=sz3_generate_prediction_positions(n0,n1,n2,n12,max_level,left_indices0,right_indices0,mid_indices0,compressed_indices0,
            left_indices1,right_indices1,mid_indices1,compressed_indices1,left_indices2,right_indices2,mid_indices2,compressed_indices2,prediction_level_offsets,dimension_order,
            reference_eb_scale<1?&reference_mask:nullptr);
    vector<uint32_t> corner_positions=sz3_generate_corner_positions(n0,n1,n2);
    vector<sz3_prediction_stencil<T>> prediction_stencils;
    if (use_position_bias)
    {
        prediction_stencils=sz3_generate_prediction_stencils(n0,n1,n2,n12,max_level,left_indices0,right_indices0,mid_indices0,compressed_indices0,left_weight0,right_weight0,cubic_stencils0,
            left_indices1,right_indices1,mid_indices1,compressed_indices1,left_weight1,right_weight1,cubic_stencils1,left_indices2,right_indices2,mid_indices2,compressed_indices2,left_weight2,right_weight2,cubic_stencils2,dimension_order,use_cubic);
        if (corner_positions.size()!=seed_num or prediction_stencils.size()!=prediction_positions.size())
            throw runtime_error("SZ3位置bias的预测位置数量错误");
        for (uint32_t i=0;i<prediction_positions.size();i++)
            if (prediction_stencils[i].target!=prediction_positions[i])
                throw runtime_error("SZ3位置bias的预测顺序不一致");
    }
    T reference_abs_eb=abs_eb*reference_eb_scale;
    T reference_neg_2_eb=1/(2*reference_abs_eb);
    if (huffman_tree_num==0)
        huffman_tree_num=min<uint32_t>(5,prediction_level_offsets.size());
    if ((use_surface_first or use_normal_boundary or use_index_order) and huffman_tree_num!=2)
        throw runtime_error("SZ3新预测树当前固定使用角点和非角点两棵Huffman tree");
    if (huffman_tree_num<1 or huffman_tree_num>5)
        throw runtime_error("SZ3 Huffman tree数量必须为1、2、3、4或5");
    vector<uint32_t> residual_group_offsets;
    if (huffman_tree_num>=3)
    {
        if (prediction_level_offsets.size()<huffman_tree_num)
            throw runtime_error("SZ3 Huffman tree数量超过当前element的角点和预测层总数");
        residual_group_offsets.push_back(0);
        for (uint32_t level=1;level<huffman_tree_num-1;level++)
            residual_group_offsets.push_back(prediction_level_offsets[level]);
        residual_group_offsets.push_back(residual_num);
    }
    double start_time=omp_get_wtime();
    vector<sz3_prediction_plan<T>> sampled_plans;
    if (sample_configuration)
    {
        if (use_surface_first or use_normal_boundary or use_index_order or use_index_midpoint or reference_eb_scale<1 or use_position_bias)
            throw runtime_error("SZ3 partition采样当前只支持标准GLL分层预测");
        auto all_orders=sz3_all_dimension_orders();
        for (uint32_t cubic_id=0;cubic_id<2;cubic_id++)
            for (auto &order:all_orders)
            {
                sz3_prediction_plan<T> plan;
                plan.dimension_order=order;
                plan.use_cubic=cubic_id==1;
                plan.positions=sz3_generate_prediction_positions(n0,n1,n2,n12,max_level,left_indices0,right_indices0,mid_indices0,compressed_indices0,
                    left_indices1,right_indices1,mid_indices1,compressed_indices1,left_indices2,right_indices2,mid_indices2,compressed_indices2,plan.level_offsets,order);
                plan.stencils=sz3_generate_prediction_stencils(n0,n1,n2,n12,max_level,left_indices0,right_indices0,mid_indices0,compressed_indices0,left_weight0,right_weight0,cubic_stencils0,
                    left_indices1,right_indices1,mid_indices1,compressed_indices1,left_weight1,right_weight1,cubic_stencils1,left_indices2,right_indices2,mid_indices2,compressed_indices2,left_weight2,right_weight2,cubic_stencils2,order,plan.use_cubic);
                if (plan.positions.size()!=plan.stencils.size() or plan.level_offsets!=prediction_level_offsets)
                    throw runtime_error("SZ3采样配置的预测层结构不一致");
                for (uint32_t i=0;i<plan.positions.size();i++)
                    if (plan.positions[i]!=plan.stencils[i].target)
                        throw runtime_error("SZ3采样配置的预测顺序不一致");
                sampled_plans.push_back(move(plan));
            }
    }
    if constexpr (avx2_rounding!=sz3_avx2_rounding::disabled)
    {
#if !defined(__AVX2__)
        throw runtime_error("当前编译未启用AVX2");
#endif
        if (!sample_configuration)
            throw runtime_error("SZ3 AVX2实验当前只支持采样优化版");
    }
    double sampling_time=0;
    vector<unsigned char> selected_plan_ids(blocks.size(),0);
    if (sample_configuration)
    {
        double sampling_start_time=start_time;
        #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
        for (uint64_t block_id=0;block_id<blocks.size();block_id++)
        {
            auto &block=blocks[block_id];
            uint64_t element_num=block.end_index-block.begin_index;
            uint64_t sample_num=min<uint64_t>(element_num,min<uint64_t>(1024,max<uint64_t>(128,element_num/4)));
            vector<T> original_sample(sample_num*n012);
            for (uint64_t sample_id=0;sample_id<sample_num;sample_id++)
            {
                uint64_t relative_element_id=((2*sample_id+1)*element_num)/(2*sample_num);
                memcpy(original_sample.data()+sample_id*n012,data.data()+(block.begin_index+relative_element_id)*n012,n012*sizeof(T));
            }
            uint32_t best_plan_id=0;
            uint64_t best_size=sz3_sample_compressed_size(original_sample,n012,corner_positions,sampled_plans[0],residual_group_offsets,huffman_tree_num,abs_eb,0);
            for (uint32_t plan_id=1;plan_id<6;plan_id++)
            {
                uint64_t candidate_size=sz3_sample_compressed_size(original_sample,n012,corner_positions,sampled_plans[plan_id],residual_group_offsets,huffman_tree_num,abs_eb,plan_id);
                if (candidate_size<best_size)
                    best_size=candidate_size,best_plan_id=plan_id;
            }
            uint32_t cubic_plan_id=best_plan_id+6;
            uint64_t cubic_size=sz3_sample_compressed_size(original_sample,n012,corner_positions,sampled_plans[cubic_plan_id],residual_group_offsets,huffman_tree_num,abs_eb,cubic_plan_id);
            if (cubic_size<best_size)
                best_plan_id=cubic_plan_id;
            selected_plan_ids[block_id]=best_plan_id;
        }
        sampling_time=omp_get_wtime()-sampling_start_time;
        vector<uint64_t> configuration_counts(sampled_plans.size(),0);
        for (auto plan_id:selected_plan_ids)
            configuration_counts[plan_id]++;
        auto all_orders=sz3_all_dimension_orders();
        printf("Sampled configurations:");
        for (uint32_t plan_id=0;plan_id<configuration_counts.size();plan_id++)
            if (configuration_counts[plan_id]>0)
                printf(" %c%c%c-%s=%lu","XYZ"[all_orders[plan_id%6][0]],"XYZ"[all_orders[plan_id%6][1]],"XYZ"[all_orders[plan_id%6][2]],plan_id<6?"linear":"cubic",configuration_counts[plan_id]);
        printf("\n");
    }
    int64_t *corner_qb=huffman_tree_num==1?nullptr:new int64_t[seed_num*params.element_num_to_process];
    int64_t *no_corner_err_qb=huffman_tree_num==2?new int64_t[params.element_num_to_process*residual_num]:nullptr;
    vector<int64_t*> residual_group_qb;
    for (uint32_t group_id=0;group_id+1<residual_group_offsets.size();group_id++)
        residual_group_qb.push_back(new int64_t[params.element_num_to_process*(residual_group_offsets[group_id+1]-residual_group_offsets[group_id])]);
    vector<uint32_t> residual_group_ids(residual_num),residual_group_local_ids(residual_num);
    for (uint32_t group_id=0;group_id<residual_group_qb.size();group_id++)
        for (uint32_t residual_id=residual_group_offsets[group_id];residual_id<residual_group_offsets[group_id+1];residual_id++)
            residual_group_ids[residual_id]=group_id,residual_group_local_ids[residual_id]=residual_id-residual_group_offsets[group_id];
    double prediction_quantization_time=0,huffman_time=0,zstd_time=0;
    uint64_t avx2_value_num=0,avx2_fallback_value_num=0,avx2_scalar_tail_value_num=0;
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus) reduction(+:prediction_quantization_time,huffman_time,zstd_time,avx2_value_num,avx2_fallback_value_num,avx2_scalar_tail_value_num)
    for (uint64_t block_id=0;block_id<blocks.size();block_id++) //这个block_id是假的，指的是在需要处理的blocks中的编号，而不是header文件中的编号
    {
        auto &block=blocks[block_id];
        uint64_t begin_index=block.begin_index;
        uint64_t end_index=block.end_index;
        uint64_t element_num=end_index-begin_index;
        vector<sz3_prediction_stencil<T>> *block_prediction_stencils=sample_configuration?&sampled_plans[selected_plan_ids[block_id]].stencils:nullptr;
        vector<int64_t> block_qb(huffman_tree_num==1?n012*element_num:0);
        int64_t *block_corner_qb=huffman_tree_num==1?nullptr:corner_qb+seed_num*begin_index;
        int64_t *block_no_corner_err_qb=huffman_tree_num==2?no_corner_err_qb+begin_index*residual_num:nullptr;
        vector<int64_t*> block_residual_group_qb;
        for (uint32_t group_id=0;group_id<residual_group_qb.size();group_id++)
            block_residual_group_qb.push_back(residual_group_qb[group_id]+begin_index*(residual_group_offsets[group_id+1]-residual_group_offsets[group_id]));
        double prediction_quantization_start_time=omp_get_wtime();
        vector<int32_t> position_bias(n012,0);
        T bias_unit=abs_eb/8;
        if (use_position_bias)
        {
            uint64_t sample_num=min<uint64_t>(element_num,min<uint64_t>(1024,max<uint64_t>(128,element_num/4)));
            vector<T> sample_data(sample_num*n012);
            for (uint64_t sample_id=0;sample_id<sample_num;sample_id++)
            {
                uint64_t relative_element_id=((2*sample_id+1)*element_num)/(2*sample_num);
                memcpy(sample_data.data()+sample_id*n012,data.data()+(begin_index+relative_element_id)*n012,n012*sizeof(T));
            }
            auto find_density_bias=[&](auto residual)
            {
                unordered_map<int64_t,uint32_t> histogram;
                for (uint64_t sample_id=0;sample_id<sample_num;sample_id++)
                {
                    int64_t bin=llround(residual(sample_id)/bias_unit);
                    histogram[bin]++;
                }
                vector<pair<int64_t,uint32_t>> sorted_histogram(histogram.begin(),histogram.end());
                sort(sorted_histogram.begin(),sorted_histogram.end());
                uint32_t default_count=0,best_count=0,window_count=0;
                int64_t best_bias=0;
                uint32_t left=0;
                for (uint32_t right=0;right<sorted_histogram.size();right++)
                {
                    window_count+=sorted_histogram[right].second;
                    if (sorted_histogram[right].first>=-8 and sorted_histogram[right].first<=8)
                        default_count+=sorted_histogram[right].second;
                    while (sorted_histogram[right].first-sorted_histogram[left].first>16)
                        window_count-=sorted_histogram[left++].second;
                    int64_t bias=sorted_histogram[right].first-8;
                    uint64_t bias_magnitude=bias<0?uint64_t(-bias):uint64_t(bias);
                    uint64_t best_magnitude=best_bias<0?uint64_t(-best_bias):uint64_t(best_bias);
                    if (window_count>best_count or (window_count==best_count and bias_magnitude<best_magnitude))
                        best_bias=bias,best_count=window_count;
                }
                uint32_t minimum_count=max<uint32_t>(8,(sample_num+49)/50);
                uint32_t minimum_improvement=max<uint32_t>(2,(sample_num+99)/100);
                if (best_count<minimum_count or best_count<default_count+minimum_improvement)
                    best_bias=0;
                if (best_bias<numeric_limits<int32_t>::min() or best_bias>numeric_limits<int32_t>::max())
                    throw runtime_error("SZ3位置bias超出int32_t范围");
                return int32_t(best_bias);
            };
            for (auto target:corner_positions)
            {
                int32_t bias=find_density_bias([&](uint64_t sample_id){return sample_data[sample_id*n012+target];});
                position_bias[target]=bias;
                for (uint64_t sample_id=0;sample_id<sample_num;sample_id++)
                {
                    T prediction=bias*bias_unit;
                    T &truth=sample_data[sample_id*n012+target];
                    int64_t qb;
                    params.quantize_and_overwrite_err(prediction,truth,qb,abs_eb,neg_2_eb);
                }
            }
            for (auto &stencil:prediction_stencils)
            {
                int32_t bias=find_density_bias([&](uint64_t sample_id)
                {
                    T *element=sample_data.data()+sample_id*n012;
                    return element[stencil.target]-sz3_predict_from_stencil(element,stencil);
                });
                position_bias[stencil.target]=bias;
                for (uint64_t sample_id=0;sample_id<sample_num;sample_id++)
                {
                    T *element=sample_data.data()+sample_id*n012;
                    T prediction=sz3_predict_from_stencil(element,stencil)+bias*bias_unit;
                    int64_t qb;
                    params.quantize_and_overwrite_err(prediction,element[stencil.target],qb,abs_eb,neg_2_eb);
                }
            }
        }
        auto corner_qb_at=[&](uint64_t element_id,uint32_t corner_id)->int64_t&
        {
            if (huffman_tree_num==1)
                return block_qb[(element_id-begin_index)*n012+corner_id];
            return corner_qb[element_id*seed_num+corner_id];
        };
        auto residual_qb_at=[&](uint64_t element_id,uint32_t residual_id)->int64_t&
        {
            if (huffman_tree_num==1)
                return block_qb[(element_id-begin_index)*n012+seed_num+residual_id];
            if (huffman_tree_num==2)
                return no_corner_err_qb[element_id*residual_num+residual_id];
            uint32_t group_id=residual_group_ids[residual_id];
            uint32_t group_size=residual_group_offsets[group_id+1]-residual_group_offsets[group_id];
            return residual_group_qb[group_id][element_id*group_size+residual_group_local_ids[residual_id]];
        };
        // 先量化位于每个element位于角上的点，覆盖原数据
        // 之后可能要改成压缩率更高的方法
        if (use_position_bias)
            for (uint64_t i=begin_index;i<end_index;i++)
                for (uint32_t corner_id=0;corner_id<seed_num;corner_id++)
                {
                    uint32_t target=corner_positions[corner_id];
                    T prediction=position_bias[target]*bias_unit;
                    params.quantize_and_overwrite_err(prediction,data[i*n012+target],corner_qb_at(i,corner_id),abs_eb,neg_2_eb);
                }
        else if (reference_eb_scale<1)
            for (uint64_t i=begin_index;i<end_index;i++)
                for (uint32_t corner_id=0;corner_id<seed_num;corner_id++)
                {
                    uint32_t target=corner_positions[corner_id];
                    params.quantize_and_overwrite(data[i*n012+target],corner_qb_at(i,corner_id),reference_abs_eb,reference_neg_2_eb);
                }
        else
            for (uint64_t i=begin_index;i<end_index;i++)
                for (uint32_t corner_id=0;corner_id<seed_num;corner_id++)
                {
                    uint32_t target=corner_positions[corner_id];
                    params.quantize_and_overwrite(data[i*n012+target],corner_qb_at(i,corner_id),abs_eb,neg_2_eb);
                }
        // 然后量化剩余的数据
        if constexpr (avx2_rounding!=sz3_avx2_rounding::disabled)
        {
#if defined(__AVX2__)
            constexpr uint32_t lane_num=is_same_v<T,float>?8:4;
            uint64_t tile_element_num=max<uint64_t>(lane_num,(256*1024)/(n012*sizeof(T)));
            tile_element_num-=tile_element_num%lane_num;
            alignas(32) int64_t vector_qb[lane_num];
            alignas(32) T reconstructed[lane_num];
            for (uint64_t tile_begin=begin_index;tile_begin<end_index;tile_begin+=tile_element_num)
            {
                uint64_t tile_end=min(end_index,tile_begin+tile_element_num);
                for (uint32_t residual_id=0;residual_id<block_prediction_stencils->size();residual_id++)
                {
                    auto &stencil=(*block_prediction_stencils)[residual_id];
                    uint64_t i=tile_begin;
                    for (;i+lane_num<=tile_end;i+=lane_num)
                    {
                        bool used_avx2=sz3_quantize_stencil_avx2<avx2_rounding>(data.data()+i*n012,n012,stencil,abs_eb,neg_2_eb,vector_qb,reconstructed);
                        if (used_avx2)
                        {
                            for (uint32_t lane=0;lane<lane_num;lane++)
                            {
                                data[(i+lane)*n012+stencil.target]=reconstructed[lane];
                                residual_qb_at(i+lane,residual_id)=vector_qb[lane];
                            }
                            avx2_value_num+=lane_num;
                        }
                        else
                        {
                            for (uint32_t lane=0;lane<lane_num;lane++)
                            {
                                T *element=data.data()+(i+lane)*n012;
                                T prediction=sz3_predict_from_stencil(element,stencil);
                                sz3_quantize_scalar<avx2_rounding>(prediction,element[stencil.target],residual_qb_at(i+lane,residual_id),abs_eb,neg_2_eb);
                            }
                            avx2_fallback_value_num+=lane_num;
                        }
                    }
                    for (;i<tile_end;i++)
                    {
                        T *element=data.data()+i*n012;
                        T prediction=sz3_predict_from_stencil(element,stencil);
                        sz3_quantize_scalar<avx2_rounding>(prediction,element[stencil.target],residual_qb_at(i,residual_id),abs_eb,neg_2_eb);
                        avx2_scalar_tail_value_num++;
                    }
                }
            }
#endif
        }
        else
            for (uint64_t i=begin_index;i<end_index;i++)
            {
                T *element=data.data()+i*n012;
                uint32_t residual_id=0;
                auto quantize_prediction=[&](T pred_value,uint32_t target)
                {
                    int64_t &qb=residual_qb_at(i,residual_id++);
                    if (use_position_bias)
                        pred_value+=position_bias[target]*bias_unit;
                    if (reference_eb_scale<1 and reference_mask[target])
                        params.quantize_and_overwrite_err(pred_value,element[target],qb,reference_abs_eb,reference_neg_2_eb);
                    else
                        params.quantize_and_overwrite_err(pred_value,element[target],qb,abs_eb,neg_2_eb);
                };
                if (sample_configuration)
                {
                    for (auto &stencil:*block_prediction_stencils)
                        quantize_prediction(sz3_predict_from_stencil(element,stencil),stencil.target);
                }
                else if (use_index_order)
                {
                    for (auto &stencil:index_order_stencils)
                    {
                        T pred_value=sz3_index_order_predict(element,stencil);
                        quantize_prediction(pred_value,stencil.target);
                    }
                }
                else if (use_normal_boundary)
                {
                    for (auto &stencil:normal_boundary_stencils)
                    {
                        T pred_value=sz3_normal_boundary_predict(element,stencil);
                        quantize_prediction(pred_value,stencil.target);
                    }
                }
                else if (use_surface_first)
                {
                    for (auto &stencil:surface_first_stencils)
                    {
                        T pred_value=sz3_surface_first_predict(element,stencil);
                        quantize_prediction(pred_value,stencil.target);
                    }
                }
                else
                {
                    uint32_t level0=0,level1=0,level2=0;
                    for (uint32_t level=0;level<max_level;level++)
                    {
                        for (auto dimension:dimension_order)
                        {
                            if (dimension==0 and n0>1 and level0<left_indices0.size())
                            {
                                for (uint32_t i0=0;i0<left_indices0[level0].size();i0++)
                                    for (auto i1:compressed_indices1[level1])
                                        for (auto i2:compressed_indices2[level2])
                                        {
                                            uint32_t pos_base=i1*n2+i2;
                                            uint32_t mid_pos=pos_base+mid_indices0[level0][i0]*n12;
                                            T pred_value=0;
                                            if (use_cubic and cubic_stencils0[level0][i0].point_num>0)
                                                for (uint32_t j=0;j<cubic_stencils0[level0][i0].point_num;j++)
                                                    pred_value+=cubic_stencils0[level0][i0].weights[j]*element[pos_base+cubic_stencils0[level0][i0].indices[j]*n12];
                                            else
                                                pred_value=left_weight0[level0][i0]*element[pos_base+left_indices0[level0][i0]*n12]+
                                                    right_weight0[level0][i0]*element[pos_base+right_indices0[level0][i0]*n12];
                                            quantize_prediction(pred_value,mid_pos);
                                        }
                                level0++;
                            }
                            else if (dimension==1 and n1>1 and level1<left_indices1.size())
                            {
                                for (auto i0:compressed_indices0[level0])
                                    for (uint32_t i1=0;i1<left_indices1[level1].size();i1++)
                                        for (auto i2:compressed_indices2[level2])
                                        {
                                            uint32_t pos_base=i0*n12+i2;
                                            uint32_t mid_pos=pos_base+mid_indices1[level1][i1]*n2;
                                            T pred_value=0;
                                            if (use_cubic and cubic_stencils1[level1][i1].point_num>0)
                                                for (uint32_t j=0;j<cubic_stencils1[level1][i1].point_num;j++)
                                                    pred_value+=cubic_stencils1[level1][i1].weights[j]*element[pos_base+cubic_stencils1[level1][i1].indices[j]*n2];
                                            else
                                                pred_value=left_weight1[level1][i1]*element[pos_base+left_indices1[level1][i1]*n2]+
                                                    right_weight1[level1][i1]*element[pos_base+right_indices1[level1][i1]*n2];
                                            quantize_prediction(pred_value,mid_pos);
                                        }
                                level1++;
                            }
                            else if (dimension==2 and n2>1 and level2<left_indices2.size())
                            {
                                for (auto i0:compressed_indices0[level0])
                                    for (auto i1:compressed_indices1[level1])
                                        for (uint32_t i2=0;i2<left_indices2[level2].size();i2++)
                                        {
                                            uint32_t pos_base=i0*n12+i1*n2;
                                            uint32_t mid_pos=pos_base+mid_indices2[level2][i2];
                                            T pred_value=0;
                                            if (use_cubic and cubic_stencils2[level2][i2].point_num>0)
                                                for (uint32_t j=0;j<cubic_stencils2[level2][i2].point_num;j++)
                                                    pred_value+=cubic_stencils2[level2][i2].weights[j]*element[pos_base+cubic_stencils2[level2][i2].indices[j]];
                                            else
                                                pred_value=left_weight2[level2][i2]*element[pos_base+left_indices2[level2][i2]]+
                                                    right_weight2[level2][i2]*element[pos_base+right_indices2[level2][i2]];
                                            quantize_prediction(pred_value,mid_pos);
                                        }
                                level2++;
                            }
                        }
                    }
                }
                if (residual_id!=residual_num)
                    throw runtime_error("SZ3预测产生的QB数量错误");
            }
        prediction_quantization_time+=omp_get_wtime()-prediction_quantization_start_time;
        // 把两种qb无损压缩到一起
        uint64_t configuration_bytes=sample_configuration?1:0;
        uint64_t position_bias_bytes=use_position_bias?n012*sizeof(int32_t):0;
        vector<unsigned char> huffman_bs_vec(n012*element_num*sizeof(int64_t)+position_bias_bytes+configuration_bytes);
        unsigned char *huffman_bs=huffman_bs_vec.data();
        if (sample_configuration)
            *huffman_bs++=selected_plan_ids[block_id];
        if (use_position_bias)
        {
            memcpy(huffman_bs,position_bias.data(),position_bias_bytes);
            huffman_bs+=position_bias_bytes;
        }
        double huffman_start_time=omp_get_wtime();
        if (huffman_tree_num==1)
        {
            mrcz::huffman_encoder_class<int64_t> encoder;
            encoder.preprocess_encode(block_qb.data(),block_qb.size());
            encoder.save(huffman_bs);
            encoder.encode(block_qb.data(),block_qb.size(),huffman_bs);
        }
        else
        {
            mrcz::huffman_encoder_class<int64_t> corner_encoder;
            corner_encoder.preprocess_encode(block_corner_qb,seed_num*element_num);
            corner_encoder.save(huffman_bs);
            corner_encoder.encode(block_corner_qb,seed_num*element_num,huffman_bs);
            if (huffman_tree_num==2)
            {
                mrcz::huffman_encoder_class<int64_t> err_encoder;
                err_encoder.preprocess_encode(block_no_corner_err_qb,residual_num*element_num);
                err_encoder.save(huffman_bs);
                err_encoder.encode(block_no_corner_err_qb,residual_num*element_num,huffman_bs);
            }
            else
            {
                for (uint32_t group_id=0;group_id<block_residual_group_qb.size();group_id++)
                {
                    uint64_t group_value_num=(residual_group_offsets[group_id+1]-residual_group_offsets[group_id])*element_num;
                    mrcz::huffman_encoder_class<int64_t> encoder;
                    encoder.preprocess_encode(block_residual_group_qb[group_id],group_value_num);
                    encoder.save(huffman_bs);
                    encoder.encode(block_residual_group_qb[group_id],group_value_num,huffman_bs);
                }
            }
        }
        huffman_time+=omp_get_wtime()-huffman_start_time;
        uint64_t huffman_bs_len=huffman_bs-huffman_bs_vec.data();
        double zstd_start_time=omp_get_wtime();
        uint64_t max_zstd_bs_len=ZSTD_compressBound(huffman_bs_len);
        block.zstd_bs.resize(max_zstd_bs_len);
        uint64_t zstd_bs_len=ZSTD_compress(block.zstd_bs.data(),max_zstd_bs_len,huffman_bs_vec.data(),huffman_bs_len,3);
        block.zstd_bs.resize(zstd_bs_len);
        zstd_time+=omp_get_wtime()-zstd_start_time;
    }
    double end_time=omp_get_wtime();
    if (sample_configuration)
        printf("Sampling time= %f seconds\n",sampling_time);
    if constexpr (avx2_rounding!=sz3_avx2_rounding::disabled)
        printf("AVX2 residual values= %lu, int64 fallback values= %lu, scalar tail values= %lu\n",avx2_value_num,avx2_fallback_value_num,avx2_scalar_tail_value_num);
    double measured_stage_time=prediction_quantization_time+huffman_time+zstd_time;
    printf("Compression accumulated partition times: prediction+quantization+QB placement= %f, Huffman= %f, Zstd= %f seconds\n",prediction_quantization_time,huffman_time,zstd_time);
    if (measured_stage_time>0)
        printf("Compression stage shares: prediction+quantization+QB placement= %.2f%%, Huffman= %.2f%%, Zstd= %.2f%%\n",100*prediction_quantization_time/measured_stage_time,100*huffman_time/measured_stage_time,100*zstd_time/measured_stage_time);
    printf("Compression time= %f seconds\n",end_time-start_time);
    printf("\033[91mField= %s, compression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),params.element_num_to_process*n012*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
    delete[] corner_qb;
    delete[] no_corner_err_qb;
    for (auto group_qb:residual_group_qb)
        delete[] group_qb;
}

template<typename T>
void sz3_3d_gll_compress(mrcz::param_class<T> &params,string field)
{
    sz3_3d_gll_compress_impl(params,field,1);
}

template<typename T>
void fdwz_compress(mrcz::param_class<T> &params,string field)
{
    sz3_3d_gll_compress_impl(params,field,0,{0,1,2},false,false,false,false,false,T(1),false,true);
}

template<typename T>
void fdwz_avx2_even_compress(mrcz::param_class<T> &params,string field)
{
    sz3_3d_gll_compress_impl<T,sz3_avx2_rounding::ties_to_even>(params,field,0,{0,1,2},false,false,false,false,false,T(1),false,true);
}

template<typename T>
void fdwz_avx2_away_compress(mrcz::param_class<T> &params,string field)
{
    sz3_3d_gll_compress_impl<T,sz3_avx2_rounding::ties_away_from_zero>(params,field,0,{0,1,2},false,false,false,false,false,T(1),false,true);
}

template<typename T>
void sz3_2_compress(mrcz::param_class<T> &params,string field)
{
    sz3_3d_gll_compress_impl(params,field,2);
}

template<typename T>
void sz3_3d_gll_experimental_compress(mrcz::param_class<T> &params,string field)
{
    sz3_3d_gll_compress_impl(params,field,5);
}

template<typename T>
void sz3_3d_gll_order_compress(mrcz::param_class<T> &params,string field)
{
    sz3_3d_gll_compress_impl(params,field,2,sz3_dimension_order_from_method(params.method));
}

template<typename T>
void sz3_3d_gll_cubic_compress(mrcz::param_class<T> &params,string field)
{
    sz3_3d_gll_compress_impl(params,field,2,{0,1,2},true);
}

template<typename T>
void sz3_3d_gll_surface_first_compress(mrcz::param_class<T> &params,string field)
{
    sz3_3d_gll_compress_impl(params,field,2,{0,1,2},false,true);
}

template<typename T>
void sz3_3d_gll_normal_boundary_compress(mrcz::param_class<T> &params,string field)
{
    sz3_3d_gll_compress_impl(params,field,2,{0,1,2},false,false,true);
}

template<typename T>
void sz3_3d_gll_index_order_compress(mrcz::param_class<T> &params,string field)
{
    sz3_3d_gll_compress_impl(params,field,2,{0,1,2},false,false,false,true);
}

template<typename T>
void sz3_3d_gll_index_midpoint_compress(mrcz::param_class<T> &params,string field)
{
    sz3_3d_gll_compress_impl(params,field,2,{0,1,2},false,false,false,false,true);
}

template<typename T>
void sz3_3d_gll_reference_eb_compress(mrcz::param_class<T> &params,string field)
{
    sz3_3d_gll_compress_impl(params,field,2,{0,1,2},false,false,false,false,false,T(sz3_reference_eb_scale_from_method(params.method)));
}

template<typename T>
void sz3_3d_gll_sampled_configuration_compress(mrcz::param_class<T> &params,string field)
{
    auto configuration=sz3_sampled_configuration_from_method(params.method);
    sz3_3d_gll_compress_impl(params,field,configuration.huffman_tree_num,configuration.dimension_order,configuration.use_cubic);
}

template<typename T>
void sz3_3d_gll_position_bias_compress(mrcz::param_class<T> &params,string field)
{
    auto configuration=sz3_position_bias_configuration_from_method(params.method);
    sz3_3d_gll_compress_impl(params,field,configuration.huffman_tree_num,configuration.dimension_order,configuration.use_cubic,false,false,false,false,T(1),true);
}

namespace mrcz
{
    template<typename T>
    inline void param_class<T>::quantize_and_overwrite_corner_numbers(int64_t *corner_qb,uint64_t begin_index,uint64_t end_index,uint32_t n0,uint32_t n1,uint32_t n2,
        uint32_t n12,uint32_t n012,T abs_eb,T neg_2_eb)
    {
        if ((n0>1) and (n1>1) and (n2>1))
            for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
            {
                T *element=data.data()+element_id*n012;
                int64_t *element_corner_qb=corner_qb+8*element_id;
                quantize_and_overwrite(element[0],element_corner_qb[0],abs_eb,neg_2_eb);
                quantize_and_overwrite(element[n2-1],element_corner_qb[1],abs_eb,neg_2_eb);
                quantize_and_overwrite(element[n12-n2],element_corner_qb[2],abs_eb,neg_2_eb);
                quantize_and_overwrite(element[n12-1],element_corner_qb[3],abs_eb,neg_2_eb);
                quantize_and_overwrite(element[n012-n12],element_corner_qb[4],abs_eb,neg_2_eb);
                quantize_and_overwrite(element[n012-n12+n2-1],element_corner_qb[5],abs_eb,neg_2_eb);
                quantize_and_overwrite(element[n012-n2],element_corner_qb[6],abs_eb,neg_2_eb);
                quantize_and_overwrite(element[n012-1],element_corner_qb[7],abs_eb,neg_2_eb);
            }
        else if ((n0>1) and (n1>1) and (n2==1))
            for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
            {
                T *element=data.data()+element_id*n012;
                int64_t *element_corner_qb=corner_qb+4*element_id;
                quantize_and_overwrite(element[0],element_corner_qb[0],abs_eb,neg_2_eb);
                quantize_and_overwrite(element[n12-n2],element_corner_qb[1],abs_eb,neg_2_eb);
                quantize_and_overwrite(element[n012-n12],element_corner_qb[2],abs_eb,neg_2_eb);
                quantize_and_overwrite(element[n012-n2],element_corner_qb[3],abs_eb,neg_2_eb);
            }
        else if ((n0>1) and (n1==1) and (n2>1))
            for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
            {
                T *element=data.data()+element_id*n012;
                int64_t *element_corner_qb=corner_qb+4*element_id;
                quantize_and_overwrite(element[0],element_corner_qb[0],abs_eb,neg_2_eb);
                quantize_and_overwrite(element[n2-1],element_corner_qb[1],abs_eb,neg_2_eb);
                quantize_and_overwrite(element[n012-n12],element_corner_qb[2],abs_eb,neg_2_eb);
                quantize_and_overwrite(element[n012-n12+n2-1],element_corner_qb[3],abs_eb,neg_2_eb);
            }
        else if ((n0==1) and (n1>1) and (n2>1))
            for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
            {
                T *element=data.data()+element_id*n012;
                int64_t *element_corner_qb=corner_qb+4*element_id;
                quantize_and_overwrite(element[0],element_corner_qb[0],abs_eb,neg_2_eb);
                quantize_and_overwrite(element[n2-1],element_corner_qb[1],abs_eb,neg_2_eb);
                quantize_and_overwrite(element[n12-n2],element_corner_qb[2],abs_eb,neg_2_eb);
                quantize_and_overwrite(element[n12-1],element_corner_qb[3],abs_eb,neg_2_eb);
            }
        else
            throw runtime_error("Error: sz3_2_compress only supports 3D or 2D elements.\n");
    }

    template<typename T>
    inline void param_class<T>::convert_err_eb_to_no_corner_err_eb(int64_t *element_no_corner_err_qb,int64_t *element_err_qb,uint32_t n0,uint32_t n1,uint32_t n2,uint32_t n12,uint32_t n012)
    {
        if ((n0>1) and (n1>1) and (n2>1))
        {
            memcpy(element_no_corner_err_qb,element_err_qb+1,sizeof(int64_t)*(n2-2));
            memcpy(element_no_corner_err_qb+n2-2,element_err_qb+n2,sizeof(int64_t)*(n12-2*n2));
            memcpy(element_no_corner_err_qb+n12-n2-2,element_err_qb+n12-n2+1,sizeof(int64_t)*(n2-2));
            memcpy(element_no_corner_err_qb+n12-4,element_err_qb+n12,sizeof(int64_t)*(n012-2*n12));
            memcpy(element_no_corner_err_qb+n012-n12-4,element_err_qb+n012-n12+1,sizeof(int64_t)*(n2-2));
            memcpy(element_no_corner_err_qb+n012-n12+n2-6,element_err_qb+n012-n12+n2,sizeof(int64_t)*(n12-2*n2));
            memcpy(element_no_corner_err_qb+n012-n2-6,element_err_qb+n012-n2+1,sizeof(int64_t)*(n2-2));
        }
        else if ((n0>1) and (n1>1) and (n2==1))
        {
            memcpy(element_no_corner_err_qb,element_err_qb+1,sizeof(int64_t)*(n12-2));
            memcpy(element_no_corner_err_qb+n12-2,element_err_qb+n12,sizeof(int64_t)*(n012-2*n12));
            memcpy(element_no_corner_err_qb+n012-n12-2,element_err_qb+n012-n12+1,sizeof(int64_t)*(n12-2));
        }
        else if ((n0>1) and (n1==1) and (n2>1))
        {
            memcpy(element_no_corner_err_qb,element_err_qb+1,sizeof(int64_t)*(n2-2));
            memcpy(element_no_corner_err_qb+n2-2,element_err_qb+n2,sizeof(int64_t)*(n012-2*n2));
            memcpy(element_no_corner_err_qb+n012-n2-2,element_err_qb+n012-n2+1,sizeof(int64_t)*(n2-2));
        }
        else if ((n0==1) and (n1>1) and (n2>1))
        {
            memcpy(element_no_corner_err_qb,element_err_qb+1,sizeof(int64_t)*(n2-2));
            memcpy(element_no_corner_err_qb+n2-2,element_err_qb+n2,sizeof(int64_t)*(n12-2*n2));
            memcpy(element_no_corner_err_qb+n12-n2-2,element_err_qb+n12-n2+1,sizeof(int64_t)*(n2-2));
        }
        else
            throw runtime_error("Error: sz3_2_compress only supports 3D or 2D elements.\n");
    }
}
