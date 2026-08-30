#pragma once

// XXX 只做预测误差诊断，不生成MRCZ bitstream。

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "mrcz.hpp"
#include "read_and_write_data.tpp"
#include "validate_params.hpp"

namespace mrcz::detail
{
    class prediction_error_stats
    {
        public:
            uint64_t count=0;
            uint64_t zero_quantized_count=0;
            long double abs_error_sum=0;
            long double squared_error_sum=0;
            vector<double> normalized_abs_errors;

            void add(double error,double data_range,double abs_eb)
            {
                double abs_error=fabs(error);
                count++;
                abs_error_sum+=abs_error;
                squared_error_sum+=error*error;
                normalized_abs_errors.push_back(abs_error/data_range);
                if (abs_error<abs_eb)
                    zero_quantized_count++;
            }

            void add_normalized(prediction_error_stats &other)
            {
                count+=other.count;
                zero_quantized_count+=other.zero_quantized_count;
                for (auto value:other.normalized_abs_errors)
                {
                    abs_error_sum+=value;
                    squared_error_sum+=value*value;
                    normalized_abs_errors.push_back(value);
                }
            }

            double percentile(double p)
            {
                if (normalized_abs_errors.empty())
                    return 0;
                vector<double> sorted_errors=normalized_abs_errors;
                sort(sorted_errors.begin(),sorted_errors.end());
                double position=p*(sorted_errors.size()-1);
                uint64_t left=position;
                uint64_t right=min<uint64_t>(left+1,sorted_errors.size()-1);
                double weight=position-left;
                return sorted_errors[left]*(1-weight)+sorted_errors[right]*weight;
            }

            void print(string name,double data_range=1,bool errors_are_normalized=false)
            {
                double mae=static_cast<double>(abs_error_sum/count);
                double rmse=sqrt(static_cast<double>(squared_error_sum/count));
                double nmae=errors_are_normalized ? mae : mae/data_range;
                double nrmse=errors_are_normalized ? rmse : rmse/data_range;
                double zero_ratio=100.0*zero_quantized_count/count;
                printf("  %s: N=%llu, MAE=%.9e, RMSE=%.9e, NMAE=%.9e, NRMSE=%.9e, zero_qb=%.3f%%\n",
                    name.c_str(),static_cast<unsigned long long>(count),mae,rmse,nmae,nrmse,zero_ratio);
                printf("    normalized |error|: p50=%.9e, p90=%.9e, p95=%.9e, p99=%.9e, max=%.9e\n",
                    percentile(0.50),percentile(0.90),percentile(0.95),percentile(0.99),
                    *max_element(normalized_abs_errors.begin(),normalized_abs_errors.end()));
            }
    };

    template<typename T>
    T linear_prediction(T x0,T x1,T value0,T value1,T target_x)
    {
        return ((x1-target_x)*value0+(target_x-x0)*value1)/(x1-x0);
    }

    inline uint64_t prediction_data_index(uint32_t i0,uint32_t i1,uint32_t i2,uint32_t n1,uint32_t n2)
    {
        return (static_cast<uint64_t>(i0)*n1+i1)*n2+i2;
    }

    inline vector<uint32_t> center_indices(uint32_t n)
    {
        if (n%2==1)
            return {n/2};
        return {n/2-1,n/2};
    }

    template<typename T>
    void analyze_prediction_field(param_class<T> &params,string field,prediction_error_stats &all_center,
        prediction_error_stats &all_corner_3,prediction_error_stats &all_corner_6)
    {
        uint32_t n0=params.header.c_element_shape[0];
        uint32_t n1=params.header.c_element_shape[1];
        uint32_t n2=params.header.c_element_shape[2];
        uint64_t element_size=params.header.element_size;
        double data_range=params.data_range_max[field]-params.data_range_min[field];
        double abs_eb=params.abs_eb[field];
        auto center0=center_indices(n0);
        auto center1=center_indices(n1);
        auto center2=center_indices(n2);
        prediction_error_stats center_stats,corner_3_stats,corner_6_stats;

        for (uint64_t element_id=0;element_id<params.element_num_to_process;element_id++)
        {
            T *element=params.data.data()+element_id*element_size;

            for (auto i0:center0)
                for (auto i1:center1)
                    for (auto i2:center2)
                    {
                        uint32_t indices[3]={i0,i1,i2};
                        T direction_predictions[3];
                        for (uint32_t dim=0;dim<3;dim++)
                        {
                            uint32_t left_indices[3]={i0,i1,i2};
                            uint32_t right_indices[3]={i0,i1,i2};
                            left_indices[dim]--;
                            right_indices[dim]++;
                            T left_value=element[prediction_data_index(left_indices[0],left_indices[1],left_indices[2],n1,n2)];
                            T right_value=element[prediction_data_index(right_indices[0],right_indices[1],right_indices[2],n1,n2)];
                            direction_predictions[dim]=linear_prediction(params.x[dim][left_indices[dim]],params.x[dim][right_indices[dim]],
                                left_value,right_value,params.x[dim][indices[dim]]);
                        }
                        T prediction=(direction_predictions[0]+direction_predictions[1]+direction_predictions[2])/3;
                        T truth=element[prediction_data_index(i0,i1,i2,n1,n2)];
                        center_stats.add(prediction-truth,data_range,abs_eb);
                    }

            for (uint32_t corner0:{0U,n0-1})
                for (uint32_t corner1:{0U,n1-1})
                    for (uint32_t corner2:{0U,n2-1})
                    {
                        uint32_t corners[3]={corner0,corner1,corner2};
                        T nearest_values[3];
                        T extrapolated_values[3];
                        for (uint32_t dim=0;dim<3;dim++)
                        {
                            uint32_t near_indices[3]={corner0,corner1,corner2};
                            uint32_t far_indices[3]={corner0,corner1,corner2};
                            if (corners[dim]==0)
                                near_indices[dim]=1,far_indices[dim]=2;
                            else
                                near_indices[dim]-=1,far_indices[dim]-=2;
                            nearest_values[dim]=element[prediction_data_index(near_indices[0],near_indices[1],near_indices[2],n1,n2)];
                            T far_value=element[prediction_data_index(far_indices[0],far_indices[1],far_indices[2],n1,n2)];
                            extrapolated_values[dim]=linear_prediction(params.x[dim][near_indices[dim]],params.x[dim][far_indices[dim]],
                                nearest_values[dim],far_value,params.x[dim][corners[dim]]);
                        }
                        T truth=element[prediction_data_index(corner0,corner1,corner2,n1,n2)];
                        T prediction_3=(nearest_values[0]+nearest_values[1]+nearest_values[2])/3;
                        T prediction_6=(extrapolated_values[0]+extrapolated_values[1]+extrapolated_values[2])/3;
                        corner_3_stats.add(prediction_3-truth,data_range,abs_eb);
                        corner_6_stats.add(prediction_6-truth,data_range,abs_eb);
                    }
        }

        printf("Prediction diagnostics for field %s, range=%.9e\n",field.c_str(),data_range);
        center_stats.print("center: 6-point, 3-direction interpolation",data_range);
        corner_3_stats.print("corner: 3-point nearest-neighbor average",data_range);
        corner_6_stats.print("corner: 6-point, 3-direction linear extrapolation",data_range);
        all_center.add_normalized(center_stats);
        all_corner_3.add_normalized(corner_3_stats);
        all_corner_6.add_normalized(corner_6_stats);
    }

    template<typename T>
    void analyze_prediction_diagnostics(param_class<T> &params)
    {
        validate_common_params(params);
        params.all_fields_to_process=params.fields_to_process.empty();
        params.read_header("compress");
        for (auto n:params.header.c_element_shape)
            if (n<4)
                throw runtime_error("预测诊断要求element的三个维度都至少包含4个点");
        params.default_init();

        prediction_error_stats all_center,all_corner_3,all_corner_6;
        for (auto &field:params.fields_to_process)
        {
            params.data.resize(params.element_num_to_process*params.header.element_size);
            ::read_data(params,field);
            auto range=minmax_element(params.data.begin(),params.data.end());
            params.data_range_min[field]=*range.first;
            params.data_range_max[field]=*range.second;
            T data_range=*range.second-*range.first;
            if (data_range==0)
            {
                printf("Prediction diagnostics skip constant field %s\n",field.c_str());
                continue;
            }
            if (params.eb_mode=="ABS")
                params.abs_eb[field]=params.input_eb,params.rel_eb[field]=params.input_eb/data_range;
            else
                params.rel_eb[field]=params.input_eb,params.abs_eb[field]=params.input_eb*data_range;
            analyze_prediction_field(params,field,all_center,all_corner_3,all_corner_6);
            params.data.clear();
            params.data.shrink_to_fit();
        }

        printf("Prediction diagnostics over all non-constant fields (normalized errors)\n");
        all_center.print("center: 6-point, 3-direction interpolation",1,true);
        all_corner_3.print("corner: 3-point nearest-neighbor average",1,true);
        all_corner_6.print("corner: 6-point, 3-direction linear extrapolation",1,true);
        printf("MRCZ prediction diagnostics complete; no bitstream was written.\n");
    }
}
