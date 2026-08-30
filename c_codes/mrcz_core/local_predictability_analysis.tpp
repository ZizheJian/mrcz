#pragma once

#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

#include "mrcz.hpp"
#include "read_and_write_data.tpp"
#include "validate_params.hpp"

namespace mrcz::detail
{
    class local_prediction_stats
    {
        public:
            uint64_t count=0;
            uint64_t zero_quantized_count=0;
            long double abs_error_sum=0;
            long double squared_error_sum=0;

            void add(double error,double abs_eb)
            {
                double abs_error=fabs(error);
                count++;
                zero_quantized_count+=llround(error/(2*abs_eb))==0;
                abs_error_sum+=abs_error;
                squared_error_sum+=error*error;
            }

            void merge(const local_prediction_stats &other)
            {
                count+=other.count;
                zero_quantized_count+=other.zero_quantized_count;
                abs_error_sum+=other.abs_error_sum;
                squared_error_sum+=other.squared_error_sum;
            }

            double mae() const
            {
                return count==0 ? 0:static_cast<double>(abs_error_sum/count);
            }

            double rmse() const
            {
                return count==0 ? 0:sqrt(static_cast<double>(squared_error_sum/count));
            }

            double zero_ratio() const
            {
                return count==0 ? 0:static_cast<double>(zero_quantized_count)/count;
            }
    };

    template<typename T>
    struct local_position_info
    {
        uint32_t indices[3];
        string node_type;
        double euclidean_distance;
        double manhattan_distance;
        double chebyshev_distance;
        bool active[3];
        uint64_t reference0[3];
        uint64_t reference1[3];
        T weight0[3];
        T weight1[3];
    };

    inline string predictability_safe_name(string value)
    {
        size_t slash=value.find_last_of("/\\");
        if (slash!=string::npos)
            value=value.substr(slash+1);
        for (char &c:value)
            if (!isalnum(static_cast<unsigned char>(c)) and c!='-' and c!='_')
                c='_';
        return value;
    }

    inline int64_t predictability_distance_key(double value)
    {
        return llround(value*1e9);
    }

    template<typename T>
    vector<local_position_info<T>> build_local_position_info(param_class<T> &params)
    {
        uint32_t n0=params.header.c_element_shape[0],n1=params.header.c_element_shape[1],n2=params.header.c_element_shape[2];
        uint32_t shape[3]={n0,n1,n2};
        uint64_t strides[3]={static_cast<uint64_t>(n1)*n2,n2,1};
        vector<local_position_info<T>> positions(params.header.element_size);

        for (uint32_t i0=0;i0<n0;i0++)
            for (uint32_t i1=0;i1<n1;i1++)
                for (uint32_t i2=0;i2<n2;i2++)
                {
                    uint32_t indices[3]={i0,i1,i2};
                    uint64_t position=(static_cast<uint64_t>(i0)*n1+i1)*n2+i2;
                    auto &info=positions[position];
                    for (int dim=0;dim<3;dim++)
                        info.indices[dim]=indices[dim];

                    int active_dimension_num=0,boundary_dimension_num=0;
                    double squared_distance=0,manhattan_distance=0,chebyshev_distance=0;
                    for (int dim=0;dim<3;dim++)
                    {
                        uint32_t n=shape[dim],index=indices[dim];
                        info.active[dim]=n>1;
                        if (n==1)
                            continue;
                        if (n<3)
                            throw runtime_error("局部可预测性诊断要求每个非退化element维度至少包含3个点");
                        active_dimension_num++;
                        boundary_dimension_num+=index==0 or index==n-1;
                        double center_distance=fabs(static_cast<double>(params.x[dim][index]));
                        squared_distance+=center_distance*center_distance;
                        manhattan_distance+=center_distance;
                        chebyshev_distance=max(chebyshev_distance,center_distance);

                        uint32_t reference0,reference1;
                        if (index==0)
                            reference0=1,reference1=2;
                        else if (index==n-1)
                            reference0=n-2,reference1=n-3;
                        else
                            reference0=index-1,reference1=index+1;
                        info.reference0[dim]=position+(static_cast<int64_t>(reference0)-index)*strides[dim];
                        info.reference1[dim]=position+(static_cast<int64_t>(reference1)-index)*strides[dim];
                        T x0=params.x[dim][reference0],x1=params.x[dim][reference1],target_x=params.x[dim][index];
                        info.weight0[dim]=(x1-target_x)/(x1-x0);
                        info.weight1[dim]=(target_x-x0)/(x1-x0);
                    }
                    info.euclidean_distance=sqrt(squared_distance);
                    info.manhattan_distance=manhattan_distance;
                    info.chebyshev_distance=chebyshev_distance;
                    if (boundary_dimension_num==active_dimension_num)
                        info.node_type="corner";
                    else if (active_dimension_num==3 and boundary_dimension_num==2)
                        info.node_type="edge";
                    else if (active_dimension_num==3 and boundary_dimension_num==1)
                        info.node_type="face";
                    else if (active_dimension_num==2 and boundary_dimension_num==1)
                        info.node_type="edge";
                    else
                        info.node_type="interior";
                }
        return positions;
    }

    inline void write_predictability_stats(ofstream &file,string scope,string group,string field,string predictor,int64_t i0,int64_t i1,int64_t i2,
        string node_type,double euclidean_distance,double manhattan_distance,double chebyshev_distance,const local_prediction_stats &stats,double data_range,double element_win_ratio=-1)
    {
        file<<scope<<','<<group<<','<<field<<','<<predictor<<','<<i0<<','<<i1<<','<<i2<<','<<node_type<<','
            <<euclidean_distance<<','<<manhattan_distance<<','<<chebyshev_distance<<','<<stats.count<<','<<stats.mae()<<','<<stats.rmse()<<','
            <<stats.mae()/data_range<<','<<stats.rmse()/data_range<<','<<stats.zero_ratio()<<','<<element_win_ratio<<'\n';
    }

    template<typename T>
    void analyze_local_predictability_field(param_class<T> &params,string field)
    {
        uint64_t element_num=params.element_num_to_process,element_size=params.header.element_size;
        if (element_num<3)
            throw runtime_error("沿element编号预测至少需要处理3个连续element");
        double data_range=params.data_range_max[field]-params.data_range_min[field];
        double abs_eb=params.abs_eb[field];
        auto positions=build_local_position_info(params);
        vector<string> predictor_names={"element"};
        for (int dim=0;dim<3;dim++)
            if (params.header.c_element_shape[dim]>1)
                predictor_names.push_back(string(1,"XYZ"[dim]));
        predictor_names.push_back("direction_mean");
        predictor_names.push_back("spatial_oracle_best");
        predictor_names.push_back("oracle_best");
        uint64_t predictor_num=predictor_names.size();
        int thread_num=max(1,min<int>(params.max_cpus,element_num));
        vector<vector<local_prediction_stats>> thread_stats(thread_num,vector<local_prediction_stats>(element_size*predictor_num));
        vector<vector<uint64_t>> thread_element_win_counts(thread_num,vector<uint64_t>(element_size));

        #pragma omp parallel for schedule(static) num_threads(thread_num)
        for (uint64_t element_id=0;element_id<element_num;element_id++)
        {
            int thread_id=omp_get_thread_num();
            T *element=params.data.data()+element_id*element_size;
            uint64_t element_reference0,element_reference1;
            T element_weight0,element_weight1;
            if (element_id==0)
                element_reference0=1,element_reference1=2,element_weight0=2,element_weight1=-1;
            else if (element_id==element_num-1)
                element_reference0=element_num-2,element_reference1=element_num-3,element_weight0=2,element_weight1=-1;
            else
                element_reference0=element_id-1,element_reference1=element_id+1,element_weight0=0.5,element_weight1=0.5;

            for (uint64_t position=0;position<element_size;position++)
            {
                auto &info=positions[position];
                T truth=element[position];
                T predictions[4];
                uint32_t direction_num=0;
                predictions[direction_num++]=element_weight0*params.data[element_reference0*element_size+position]+
                    element_weight1*params.data[element_reference1*element_size+position];
                for (int dim=0;dim<3;dim++)
                    if (info.active[dim])
                        predictions[direction_num++]=info.weight0[dim]*element[info.reference0[dim]]+info.weight1[dim]*element[info.reference1[dim]];

                T prediction_sum=0,best_prediction=predictions[0],best_spatial_prediction=predictions[1];
                for (uint32_t direction=0;direction<direction_num;direction++)
                {
                    prediction_sum+=predictions[direction];
                    if (fabs(predictions[direction]-truth)<fabs(best_prediction-truth))
                        best_prediction=predictions[direction];
                    if (direction>0 and fabs(predictions[direction]-truth)<fabs(best_spatial_prediction-truth))
                        best_spatial_prediction=predictions[direction];
                    thread_stats[thread_id][position*predictor_num+direction].add(predictions[direction]-truth,abs_eb);
                }
                thread_stats[thread_id][position*predictor_num+direction_num].add(prediction_sum/direction_num-truth,abs_eb);
                thread_stats[thread_id][position*predictor_num+direction_num+1].add(best_spatial_prediction-truth,abs_eb);
                thread_stats[thread_id][position*predictor_num+direction_num+2].add(best_prediction-truth,abs_eb);
                thread_element_win_counts[thread_id][position]+=fabs(predictions[0]-truth)<fabs(best_spatial_prediction-truth);
            }
        }

        vector<local_prediction_stats> position_stats(element_size*predictor_num);
        for (int thread_id=0;thread_id<thread_num;thread_id++)
            for (uint64_t i=0;i<position_stats.size();i++)
                position_stats[i].merge(thread_stats[thread_id][i]);
        vector<uint64_t> element_win_counts(element_size);
        for (int thread_id=0;thread_id<thread_num;thread_id++)
            for (uint64_t position=0;position<element_size;position++)
                element_win_counts[position]+=thread_element_win_counts[thread_id][position];

        vector<local_prediction_stats> overall_stats(predictor_num);
        map<string,vector<local_prediction_stats>> node_type_stats;
        map<int64_t,vector<local_prediction_stats>> euclidean_stats,manhattan_stats,chebyshev_stats;
        for (uint64_t position=0;position<element_size;position++)
        {
            auto &info=positions[position];
            auto &type_values=node_type_stats.try_emplace(info.node_type,predictor_num).first->second;
            auto &euclidean_values=euclidean_stats.try_emplace(predictability_distance_key(info.euclidean_distance),predictor_num).first->second;
            auto &manhattan_values=manhattan_stats.try_emplace(predictability_distance_key(info.manhattan_distance),predictor_num).first->second;
            auto &chebyshev_values=chebyshev_stats.try_emplace(predictability_distance_key(info.chebyshev_distance),predictor_num).first->second;
            for (uint64_t predictor=0;predictor<predictor_num;predictor++)
            {
                auto &values=position_stats[position*predictor_num+predictor];
                overall_stats[predictor].merge(values);
                type_values[predictor].merge(values);
                euclidean_values[predictor].merge(values);
                manhattan_values[predictor].merge(values);
                chebyshev_values[predictor].merge(values);
            }
        }

        string output_path="e25_"+predictability_safe_name(params.data_path)+"_"+predictability_safe_name(field)+"_predictability.csv";
        ofstream file(output_path);
        if (!file)
            throw runtime_error("无法创建局部可预测性结果文件: "+output_path);
        file<<setprecision(17);
        file<<"scope,group,field,predictor,i0,i1,i2,node_type,euclidean_distance,manhattan_distance,chebyshev_distance,N,MAE,RMSE,NMAE,NRMSE,zero_qb_ratio,element_win_ratio\n";
        for (uint64_t position=0;position<element_size;position++)
        {
            auto &info=positions[position];
            for (uint64_t predictor=0;predictor<predictor_num;predictor++)
                write_predictability_stats(file,"position","",field,predictor_names[predictor],info.indices[0],info.indices[1],info.indices[2],info.node_type,
                    info.euclidean_distance,info.manhattan_distance,info.chebyshev_distance,position_stats[position*predictor_num+predictor],data_range,
                    static_cast<double>(element_win_counts[position])/element_num);
        }
        for (uint64_t predictor=0;predictor<predictor_num;predictor++)
            write_predictability_stats(file,"overall","all",field,predictor_names[predictor],-1,-1,-1,"",-1,-1,-1,overall_stats[predictor],data_range);
        for (auto &[group,values]:node_type_stats)
            for (uint64_t predictor=0;predictor<predictor_num;predictor++)
                write_predictability_stats(file,"node_type",group,field,predictor_names[predictor],-1,-1,-1,group,-1,-1,-1,values[predictor],data_range);
        auto write_distance_groups=[&](string scope,map<int64_t,vector<local_prediction_stats>> &groups)
        {
            for (auto &[key,values]:groups)
                for (uint64_t predictor=0;predictor<predictor_num;predictor++)
                {
                    double distance=key/1e9;
                    double euclidean=scope=="euclidean_distance" ? distance:-1;
                    double manhattan=scope=="manhattan_distance" ? distance:-1;
                    double chebyshev=scope=="chebyshev_distance" ? distance:-1;
                    write_predictability_stats(file,scope,to_string(distance),field,predictor_names[predictor],-1,-1,-1,"",euclidean,manhattan,chebyshev,values[predictor],data_range);
                }
        };
        write_distance_groups("euclidean_distance",euclidean_stats);
        write_distance_groups("manhattan_distance",manhattan_stats);
        write_distance_groups("chebyshev_distance",chebyshev_stats);
        file.close();

        printf("Local predictability for field %s, selected range=[%.9e,%.9e], abs_eb=%.9e\n",field.c_str(),params.data_range_min[field],params.data_range_max[field],abs_eb);
        for (uint64_t predictor=0;predictor<predictor_num;predictor++)
            printf("  %s: N=%llu, NMAE=%.9e, NRMSE=%.9e, zero_qb=%.3f%%\n",predictor_names[predictor].c_str(),
                static_cast<unsigned long long>(overall_stats[predictor].count),overall_stats[predictor].mae()/data_range,
                overall_stats[predictor].rmse()/data_range,100*overall_stats[predictor].zero_ratio());
        uint64_t total_element_wins=0;
        for (auto count:element_win_counts)
            total_element_wins+=count;
        double spatial_oracle_nmae=overall_stats[predictor_num-2].mae()/data_range;
        double all_oracle_nmae=overall_stats[predictor_num-1].mae()/data_range;
        printf("  cross-element unique-win ratio=%.3f%%, oracle NMAE improvement over spatial-only=%.3f%%\n",
            100.0*total_element_wins/(element_num*element_size),100*(spatial_oracle_nmae-all_oracle_nmae)/spatial_oracle_nmae);
        printf("  detailed result: %s\n",output_path.c_str());
    }

    template<typename T>
    void analyze_local_predictability_loaded_header(param_class<T> &params)
    {
        params.default_init();
        params.max_cpus=get_max_cpus(static_cast<int>(min<uint64_t>(params.max_cpus,params.element_num_to_process)));
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
                printf("Local predictability skips constant field %s\n",field.c_str());
                continue;
            }
            if (params.eb_mode=="ABS")
                params.abs_eb[field]=params.input_eb,params.rel_eb[field]=params.input_eb/data_range;
            else
                params.rel_eb[field]=params.input_eb,params.abs_eb[field]=params.input_eb*data_range;
            analyze_local_predictability_field(params,field);
            params.data.clear();
            params.data.shrink_to_fit();
        }
        printf("SZ3-3D-GLL-predictability diagnostics complete; no bitstream was written.\n");
    }

    template<typename T>
    void analyze_local_predictability(param_class<T> &params)
    {
        validate_common_params(params);
        params.all_fields_to_process=params.fields_to_process.empty();
        params.read_header("compress");
        if (params.header.word_size==8 and is_same_v<T,float>)
        {
            param_class<double> double_params(params);
            analyze_local_predictability_loaded_header(double_params);
        }
        else
            analyze_local_predictability_loaded_header(params);
    }
}
