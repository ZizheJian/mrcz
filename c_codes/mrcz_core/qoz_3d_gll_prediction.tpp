#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

namespace mrcz::detail
{
    template<typename T>
    class qoz_axis_prediction_plan
    {
        public:
            class level_class
            {
                public:
                    vector<uint32_t> old_indices;
                    vector<uint32_t> left_indices;
                    vector<uint32_t> right_indices;
                    vector<uint32_t> mid_indices;
                    vector<T> left_weights;
                    vector<T> right_weights;
            };

            vector<level_class> levels;
            vector<uint32_t> final_indices;

            qoz_axis_prediction_plan(vector<T> &x)
            {
                if (x.size()==1)
                {
                    final_indices={0};
                    return;
                }

                vector<array<uint32_t,2>> intervals={{0,static_cast<uint32_t>(x.size()-1)}};
                final_indices={0,static_cast<uint32_t>(x.size()-1)};
                while (!intervals.empty())
                {
                    level_class level;
                    level.old_indices=final_indices;
                    vector<array<uint32_t,2>> next_intervals;
                    for (auto interval:intervals)
                    {
                        uint32_t left=interval[0];
                        uint32_t right=interval[1];
                        uint32_t mid=left+1;
                        T ideal_mid=(x[left]+x[right])/2;
                        for (uint32_t i=left+2;i<right;i++)
                            if (fabs(x[i]-ideal_mid)<fabs(x[mid]-ideal_mid))
                                mid=i;

                        level.left_indices.push_back(left);
                        level.right_indices.push_back(right);
                        level.mid_indices.push_back(mid);
                        level.left_weights.push_back((x[right]-x[mid])/(x[right]-x[left]));
                        level.right_weights.push_back((x[mid]-x[left])/(x[right]-x[left]));
                        if (mid-left>1)
                            next_intervals.push_back({left,mid});
                        if (right-mid>1)
                            next_intervals.push_back({mid,right});
                    }
                    final_indices.insert(final_indices.end(),level.mid_indices.begin(),level.mid_indices.end());
                    sort(final_indices.begin(),final_indices.end());
                    levels.push_back(level);
                    intervals=next_intervals;
                }
            }

            vector<uint32_t> &old_indices(uint32_t level)
            {
                if (level<levels.size())
                    return levels[level].old_indices;
                return final_indices;
            }
    };

    template<typename T>
    class qoz_3d_gll_prediction_plan
    {
        public:
            array<qoz_axis_prediction_plan<T>,3> axes;
            uint32_t level_num;

            qoz_3d_gll_prediction_plan(vector<vector<T>> &x):axes{
                qoz_axis_prediction_plan<T>(x[0]),qoz_axis_prediction_plan<T>(x[1]),qoz_axis_prediction_plan<T>(x[2])}
            {
                level_num=max({axes[0].levels.size(),axes[1].levels.size(),axes[2].levels.size()});
            }
    };

    inline uint32_t qoz_3d_gll_position(uint32_t i0,uint32_t i1,uint32_t i2,uint32_t n1,uint32_t n2)
    {
        return (i0*n1+i1)*n2+i2;
    }

    template<typename T>
    T qoz_3d_gll_direction_prediction(T *element,uint32_t indices[3],uint32_t dim,
        typename qoz_axis_prediction_plan<T>::level_class &level,uint32_t entry,uint32_t n1,uint32_t n2)
    {
        uint32_t left_indices[3]={indices[0],indices[1],indices[2]};
        uint32_t right_indices[3]={indices[0],indices[1],indices[2]};
        left_indices[dim]=level.left_indices[entry];
        right_indices[dim]=level.right_indices[entry];
        T left_value=element[qoz_3d_gll_position(left_indices[0],left_indices[1],left_indices[2],n1,n2)];
        T right_value=element[qoz_3d_gll_position(right_indices[0],right_indices[1],right_indices[2],n1,n2)];
        return level.left_weights[entry]*left_value+level.right_weights[entry]*right_value;
    }

    template<typename T,typename ApplyPrediction>
    void qoz_3d_gll_process_element(T *element,int64_t *element_err_qb,qoz_3d_gll_prediction_plan<T> &plan,
        uint32_t n1,uint32_t n2,ApplyPrediction apply_prediction)
    {
        for (uint32_t level_id=0;level_id<plan.level_num;level_id++)
        {
            auto &old0=plan.axes[0].old_indices(level_id);
            auto &old1=plan.axes[1].old_indices(level_id);
            auto &old2=plan.axes[2].old_indices(level_id);
            bool new0=level_id<plan.axes[0].levels.size();
            bool new1=level_id<plan.axes[1].levels.size();
            bool new2=level_id<plan.axes[2].levels.size();

            // One new coordinate: ordinary one-dimensional interpolation.
            if (new0)
            {
                auto &level0=plan.axes[0].levels[level_id];
                for (uint32_t e0=0;e0<level0.mid_indices.size();e0++)
                    for (auto i1:old1)
                        for (auto i2:old2)
                        {
                            uint32_t indices[3]={level0.mid_indices[e0],i1,i2};
                            uint32_t pos=qoz_3d_gll_position(indices[0],indices[1],indices[2],n1,n2);
                            T prediction=qoz_3d_gll_direction_prediction(element,indices,0,level0,e0,n1,n2);
                            apply_prediction(prediction,element[pos],element_err_qb[pos]);
                        }
            }
            if (new1)
            {
                auto &level1=plan.axes[1].levels[level_id];
                for (auto i0:old0)
                    for (uint32_t e1=0;e1<level1.mid_indices.size();e1++)
                        for (auto i2:old2)
                        {
                            uint32_t indices[3]={i0,level1.mid_indices[e1],i2};
                            uint32_t pos=qoz_3d_gll_position(indices[0],indices[1],indices[2],n1,n2);
                            T prediction=qoz_3d_gll_direction_prediction(element,indices,1,level1,e1,n1,n2);
                            apply_prediction(prediction,element[pos],element_err_qb[pos]);
                        }
            }
            if (new2)
            {
                auto &level2=plan.axes[2].levels[level_id];
                for (auto i0:old0)
                    for (auto i1:old1)
                        for (uint32_t e2=0;e2<level2.mid_indices.size();e2++)
                        {
                            uint32_t indices[3]={i0,i1,level2.mid_indices[e2]};
                            uint32_t pos=qoz_3d_gll_position(indices[0],indices[1],indices[2],n1,n2);
                            T prediction=qoz_3d_gll_direction_prediction(element,indices,2,level2,e2,n1,n2);
                            apply_prediction(prediction,element[pos],element_err_qb[pos]);
                        }
            }

            // Two new coordinates: average the two directional predictions.
            if (new0 and new1)
            {
                auto &level0=plan.axes[0].levels[level_id];
                auto &level1=plan.axes[1].levels[level_id];
                for (uint32_t e0=0;e0<level0.mid_indices.size();e0++)
                    for (uint32_t e1=0;e1<level1.mid_indices.size();e1++)
                        for (auto i2:old2)
                        {
                            uint32_t indices[3]={level0.mid_indices[e0],level1.mid_indices[e1],i2};
                            T prediction0=qoz_3d_gll_direction_prediction(element,indices,0,level0,e0,n1,n2);
                            T prediction1=qoz_3d_gll_direction_prediction(element,indices,1,level1,e1,n1,n2);
                            uint32_t pos=qoz_3d_gll_position(indices[0],indices[1],indices[2],n1,n2);
                            T prediction=(prediction0+prediction1)/2;
                            apply_prediction(prediction,element[pos],element_err_qb[pos]);
                        }
            }
            if (new0 and new2)
            {
                auto &level0=plan.axes[0].levels[level_id];
                auto &level2=plan.axes[2].levels[level_id];
                for (uint32_t e0=0;e0<level0.mid_indices.size();e0++)
                    for (auto i1:old1)
                        for (uint32_t e2=0;e2<level2.mid_indices.size();e2++)
                        {
                            uint32_t indices[3]={level0.mid_indices[e0],i1,level2.mid_indices[e2]};
                            T prediction0=qoz_3d_gll_direction_prediction(element,indices,0,level0,e0,n1,n2);
                            T prediction2=qoz_3d_gll_direction_prediction(element,indices,2,level2,e2,n1,n2);
                            uint32_t pos=qoz_3d_gll_position(indices[0],indices[1],indices[2],n1,n2);
                            T prediction=(prediction0+prediction2)/2;
                            apply_prediction(prediction,element[pos],element_err_qb[pos]);
                        }
            }
            if (new1 and new2)
            {
                auto &level1=plan.axes[1].levels[level_id];
                auto &level2=plan.axes[2].levels[level_id];
                for (auto i0:old0)
                    for (uint32_t e1=0;e1<level1.mid_indices.size();e1++)
                        for (uint32_t e2=0;e2<level2.mid_indices.size();e2++)
                        {
                            uint32_t indices[3]={i0,level1.mid_indices[e1],level2.mid_indices[e2]};
                            T prediction1=qoz_3d_gll_direction_prediction(element,indices,1,level1,e1,n1,n2);
                            T prediction2=qoz_3d_gll_direction_prediction(element,indices,2,level2,e2,n1,n2);
                            uint32_t pos=qoz_3d_gll_position(indices[0],indices[1],indices[2],n1,n2);
                            T prediction=(prediction1+prediction2)/2;
                            apply_prediction(prediction,element[pos],element_err_qb[pos]);
                        }
            }

            // Three new coordinates: average all three directional predictions.
            if (new0 and new1 and new2)
            {
                auto &level0=plan.axes[0].levels[level_id];
                auto &level1=plan.axes[1].levels[level_id];
                auto &level2=plan.axes[2].levels[level_id];
                for (uint32_t e0=0;e0<level0.mid_indices.size();e0++)
                    for (uint32_t e1=0;e1<level1.mid_indices.size();e1++)
                        for (uint32_t e2=0;e2<level2.mid_indices.size();e2++)
                        {
                            uint32_t indices[3]={level0.mid_indices[e0],level1.mid_indices[e1],level2.mid_indices[e2]};
                            T prediction0=qoz_3d_gll_direction_prediction(element,indices,0,level0,e0,n1,n2);
                            T prediction1=qoz_3d_gll_direction_prediction(element,indices,1,level1,e1,n1,n2);
                            T prediction2=qoz_3d_gll_direction_prediction(element,indices,2,level2,e2,n1,n2);
                            uint32_t pos=qoz_3d_gll_position(indices[0],indices[1],indices[2],n1,n2);
                            T prediction=(prediction0+prediction1+prediction2)/3;
                            apply_prediction(prediction,element[pos],element_err_qb[pos]);
                        }
            }
        }
    }
}
