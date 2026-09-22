#pragma once

#include <array>
#include <memory>
#include "mrcz.hpp"
#include "fdwz_common.hpp"
#include "huffman.hpp"
#include "quantize.tpp"
#include <fstream>
#include <iomanip>
#include <zstd.h>

namespace mrcz::detail
{
    constexpr long double fdwz_prediction_error_bin_width=1e-3L;

    inline string fdwz_tree_group_name(uint32_t tree_id,uint32_t huffman_tree_num)
    {
        if(huffman_tree_num==1)
            return "all";
        return tree_id==0?"corner":(tree_id+1==huffman_tree_num?"remaining":"prediction-level-"+to_string(tree_id-1));
    }

    template<typename T>
    inline int32_t fdwz_prediction_error_bin(T error)
    {
        long double bin=floor(static_cast<long double>(error)/fdwz_prediction_error_bin_width+0.5L);
        if(!isfinite(bin) or bin<=numeric_limits<int32_t>::min() or bin>numeric_limits<int32_t>::max())
            return numeric_limits<int32_t>::min();
        return static_cast<int32_t>(bin);
    }

    // 复用候选方案评估所需的临时数组，避免每个partition连续评估6个方案时反复申请大块内存。
    template<typename T>
    struct fdwz_sampling_workspace
    {
        vector<T> original_sample,sample;
        vector<vector<int32_t>> tree_qb;
        vector<unsigned char> huffman_bs,zstd_bs;
    };

    // INT32_MIN作为unpredictable标记；量化后的QB无法放入int32时原值无损保存。
    template<bool check_qb_range,typename T>
    inline int32_t fdwz_quantize_or_save_unpredictable(T &value,T prediction,T abs_eb,T inv_2_eb,vector<T> &unpredictable)
    {
        T scaled=(prediction-value)*inv_2_eb;
        T rounded=round(scaled);
        if constexpr(!check_qb_range)
        {
            int32_t qb=static_cast<int32_t>(rounded);
            value=prediction-qb*2*abs_eb;
            return qb;
        }
        else
        {
            long double rounded_bound=static_cast<long double>(rounded);
            if(rounded_bound>=-static_cast<long double>(numeric_limits<int32_t>::max()) and rounded_bound<=static_cast<long double>(numeric_limits<int32_t>::max()))
            {
                int32_t qb=static_cast<int32_t>(rounded);
                value=prediction-qb*2*abs_eb;
                return qb;
            }
        }
        unpredictable.push_back(value);
        return numeric_limits<int32_t>::min();
    }

    // 根据field值域、ABS EB和全部候选stencil权重估计预测误差QB上界；无法证明int32安全时改为逐点检查。角点相对零量化，因此始终逐点检查。
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

    // 在抽样element上执行预测、量化、Huffman和Zstd，以最终字节数评价一个候选plan。
    template<typename T>
    uint64_t fdwz_sample_compressed_size(vector<T> &original_sample,uint32_t n012,vector<uint32_t> &corner_positions,fdwz_prediction_plan<T> &plan,
        vector<uint32_t> &residual_group_offsets,uint32_t huffman_tree_num,T abs_eb,unsigned char configuration,fdwz_sampling_workspace<T> &workspace)
    {
        uint64_t sample_num=original_sample.size()/n012;
        uint32_t corner_num=corner_positions.size();
        uint32_t residual_num=n012-corner_num;
        T inv_2_eb=1/(2*abs_eb);
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
                element_corner_qb[corner_id]=fdwz_quantize_or_save_unpredictable<true>(element[target],T(0),abs_eb,inv_2_eb,unpredictable);
            }
            if (huffman_tree_num==1)
            {
                int32_t *element_residual_qb=element_corner_qb+corner_num;
                for (uint32_t residual_id=0;residual_id<residual_num;residual_id++)
                {
                    auto &stencil=plan.stencils[residual_id];
                    T prediction=fdwz_predict_from_stencil(element,stencil);
                    element_residual_qb[residual_id]=fdwz_quantize_or_save_unpredictable<true>(element[stencil.target],prediction,abs_eb,inv_2_eb,unpredictable);
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
                        element_group_qb[local_id]=fdwz_quantize_or_save_unpredictable<true>(element[stencil.target],prediction,abs_eb,inv_2_eb,unpredictable);
                    }
                }
            }
        }
        // sample评分使用保守的临时空间；首字节保存当前候选plan id，供Zstd评分包含相同的metadata开销。
        auto &huffman_bs_vec=workspace.huffman_bs;
        huffman_bs_vec.resize(sample_num*n012*32+4096*huffman_tree_num+sizeof(uint64_t)+1);
        unsigned char *huffman_bs=huffman_bs_vec.data();
        *huffman_bs++=configuration|0x80;
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

    // 压缩一个partition；QB固定为int32，plan字节的bit 7标记32位QB，INT32_MIN直接标记unpredictable原值。
    template<bool measure_phase_time,bool check_residual_qb_range,bool analyze_tree_distributions,typename T>
    void fdwz_compress_block(mrcz::param_class<T> &params,uint64_t block_id,uint32_t n012,uint32_t corner_num,uint32_t residual_num,vector<uint32_t> &corner_positions,
        vector<fdwz_prediction_stencil<T>> &prediction_stencils,vector<uint32_t> &residual_group_offsets,array<uint32_t,5> &tree_value_num,uint32_t huffman_tree_num,T abs_eb,T inv_2_eb,
        unsigned char selected_plan_id,int32_t *block_qb,int32_t *block_prediction_error_bins,vector<uint64_t> &unpredictable_counts,vector<double> &prediction_quantization_time,vector<double> &huffman_time,vector<double> &zstd_time)
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
            int32_t *element_corner_prediction_error_bins=block_prediction_error_bins+relative_element_id*(huffman_tree_num==1?n012:corner_num);
            for(uint32_t corner_id=0;corner_id<corner_num;corner_id++)
            {
                uint32_t target=corner_positions[corner_id];
                if constexpr(analyze_tree_distributions)
                    element_corner_prediction_error_bins[corner_id]=fdwz_prediction_error_bin(data[i*n012+target]);
                element_corner_qb[corner_id]=fdwz_quantize_or_save_unpredictable<true>(data[i*n012+target],T(0),abs_eb,inv_2_eb,unpredictable);
            }
        }
        if(huffman_tree_num==1)
        {
            for(uint64_t i=begin_index;i<end_index;i++)
            {
                T *element=data.data()+i*n012;
                uint64_t relative_element_id=i-begin_index;
                int32_t *element_residual_qb=block_qb+relative_element_id*n012+corner_num;
                int32_t *element_residual_prediction_error_bins=block_prediction_error_bins+relative_element_id*n012+corner_num;
                for(uint32_t residual_id=0;residual_id<residual_num;residual_id++)
                {
                    auto &stencil=prediction_stencils[residual_id];
                    T prediction=fdwz_predict_from_stencil(element,stencil);
                    if constexpr(analyze_tree_distributions)
                        element_residual_prediction_error_bins[residual_id]=fdwz_prediction_error_bin(element[stencil.target]-prediction);
                    element_residual_qb[residual_id]=fdwz_quantize_or_save_unpredictable<check_residual_qb_range>(element[stencil.target],prediction,abs_eb,inv_2_eb,unpredictable);
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
                    int32_t *batch_group_prediction_error_bins=block_prediction_error_bins+tree_offsets[group_id+1]+batch_relative_element_id*group_size;
                    for(uint32_t local_id=0;local_id<group_size;local_id++)
                    {
                        auto &stencil=prediction_stencils[group_begin+local_id];
                        #pragma GCC unroll 4
                        for(uint64_t lane=0;lane<batch_size;lane++)
                        {
                            T *batch_element=data.data()+(batch_begin+lane)*n012;
                            T prediction=fdwz_predict_from_stencil(batch_element,stencil);
                            if constexpr(analyze_tree_distributions)
                                batch_group_prediction_error_bins[lane*group_size+local_id]=fdwz_prediction_error_bin(batch_element[stencil.target]-prediction);
                            batch_group_qb[lane*group_size+local_id]=fdwz_quantize_or_save_unpredictable<check_residual_qb_range>(batch_element[stencil.target],prediction,abs_eb,inv_2_eb,unpredictable);
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
                    int32_t *element_group_prediction_error_bins=block_prediction_error_bins+tree_offsets[group_id+1]+relative_element_id*group_size;
                    for(uint32_t local_id=0;local_id<group_size;local_id++)
                    {
                        auto &stencil=prediction_stencils[group_begin+local_id];
                        T prediction=fdwz_predict_from_stencil(element,stencil);
                        if constexpr(analyze_tree_distributions)
                            element_group_prediction_error_bins[local_id]=fdwz_prediction_error_bin(element[stencil.target]-prediction);
                        element_group_qb[local_id]=fdwz_quantize_or_save_unpredictable<check_residual_qb_range>(element[stencil.target],prediction,abs_eb,inv_2_eb,unpredictable);
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
        *huffman_bs++=selected_plan_id|0x80;
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

    // 导出正式量化阶段中每棵Huffman tree实际接收的QB频数；统计在压缩计时结束后执行，不改变bitstream。
    template<typename T>
    void fdwz_write_tree_qb_distribution(mrcz::param_class<T> &params,string field,uint32_t huffman_tree_num,array<uint32_t,5> &tree_value_num,vector<unique_ptr<int32_t[]>> &partition_qb)
    {
        vector<unordered_map<int32_t,uint64_t>> distributions(huffman_tree_num);
        for(uint64_t block_id=0;block_id<params.partitions.size();block_id++)
        {
            uint64_t element_num=params.partitions[block_id].element_num;
            uint64_t offset=0;
            for(uint32_t tree_id=0;tree_id<huffman_tree_num;tree_id++)
            {
                uint64_t tree_qb_num=tree_value_num[tree_id]*element_num;
                int32_t *tree_qb=partition_qb[block_id].get()+offset;
                auto &distribution=distributions[tree_id];
                for(uint64_t qb_id=0;qb_id<tree_qb_num;qb_id++)
                    distribution[tree_qb[qb_id]]++;
                offset+=tree_qb_num;
            }
        }
        ifstream existing_file(params.fdwz_tree_qb_distribution_path,ios::binary|ios::ate);
        bool write_header=!existing_file.good() or existing_file.tellg()==0;
        existing_file.close();
        ofstream output(params.fdwz_tree_qb_distribution_path,ios::app);
        if(!output)
            throw runtime_error("无法写入FDWZ tree QB分布: "+params.fdwz_tree_qb_distribution_path);
        if(write_header)
            output<<"Field,REL EB,ABS EB,Tree ID,Group,QB,Count\n";
        output<<setprecision(17);
        for(uint32_t tree_id=0;tree_id<huffman_tree_num;tree_id++)
        {
            string group=fdwz_tree_group_name(tree_id,huffman_tree_num);
            vector<pair<int32_t,uint64_t>> sorted_distribution(distributions[tree_id].begin(),distributions[tree_id].end());
            sort(sorted_distribution.begin(),sorted_distribution.end(),[](auto &left,auto &right){return left.first<right.first;});
            for(auto &[qb,count]:sorted_distribution)
                output<<field<<','<<params.rel_eb[field]<<','<<params.abs_eb[field]<<','<<tree_id<<','<<group<<','<<qb<<','<<count<<'\n';
        }
        printf("FDWZ tree QB distribution written to %s\n",params.fdwz_tree_qb_distribution_path.c_str());
    }

    template<typename T>
    void fdwz_write_tree_prediction_error_distribution(mrcz::param_class<T> &params,string field,uint32_t huffman_tree_num,array<uint32_t,5> &tree_value_num,vector<unique_ptr<int32_t[]>> &partition_prediction_error_bins)
    {
        vector<int32_t> minimum_bins(huffman_tree_num,numeric_limits<int32_t>::max()),maximum_bins(huffman_tree_num,numeric_limits<int32_t>::min());
        for(uint64_t block_id=0;block_id<params.partitions.size();block_id++)
        {
            uint64_t element_num=params.partitions[block_id].element_num,offset=0;
            for(uint32_t tree_id=0;tree_id<huffman_tree_num;tree_id++)
            {
                uint64_t tree_value_count=tree_value_num[tree_id]*element_num;
                int32_t *tree_bins=partition_prediction_error_bins[block_id].get()+offset;
                for(uint64_t value_id=0;value_id<tree_value_count;value_id++)
                {
                    if(tree_bins[value_id]==numeric_limits<int32_t>::min())
                        throw runtime_error("FDWZ prediction error无法用1e-3宽度的int32 bin表示");
                    minimum_bins[tree_id]=min(minimum_bins[tree_id],tree_bins[value_id]);
                    maximum_bins[tree_id]=max(maximum_bins[tree_id],tree_bins[value_id]);
                }
                offset+=tree_value_count;
            }
        }
        vector<vector<uint64_t>> distributions(huffman_tree_num);
        for(uint32_t tree_id=0;tree_id<huffman_tree_num;tree_id++)
            distributions[tree_id].resize(static_cast<uint64_t>(maximum_bins[tree_id])-minimum_bins[tree_id]+1,0);
        for(uint64_t block_id=0;block_id<params.partitions.size();block_id++)
        {
            uint64_t element_num=params.partitions[block_id].element_num,offset=0;
            for(uint32_t tree_id=0;tree_id<huffman_tree_num;tree_id++)
            {
                uint64_t tree_value_count=tree_value_num[tree_id]*element_num;
                int32_t *tree_bins=partition_prediction_error_bins[block_id].get()+offset;
                for(uint64_t value_id=0;value_id<tree_value_count;value_id++)
                    distributions[tree_id][static_cast<int64_t>(tree_bins[value_id])-minimum_bins[tree_id]]++;
                offset+=tree_value_count;
            }
        }
        ifstream existing_file(params.fdwz_tree_prediction_error_distribution_path,ios::binary|ios::ate);
        bool write_header=!existing_file.good() or existing_file.tellg()==0;
        existing_file.close();
        ofstream output(params.fdwz_tree_prediction_error_distribution_path,ios::app);
        if(!output)
            throw runtime_error("无法写入FDWZ tree prediction-error分布: "+params.fdwz_tree_prediction_error_distribution_path);
        if(write_header)
            output<<"Field,REL EB,ABS EB,Tree ID,Group,Bin ID,Bin Width,Prediction Error,Count\n";
        output<<setprecision(17);
        for(uint32_t tree_id=0;tree_id<huffman_tree_num;tree_id++)
        {
            string group=fdwz_tree_group_name(tree_id,huffman_tree_num);
            for(uint64_t distribution_id=0;distribution_id<distributions[tree_id].size();distribution_id++)
            {
                uint64_t count=distributions[tree_id][distribution_id];
                if(count==0)
                    continue;
                int64_t bin_id=static_cast<int64_t>(minimum_bins[tree_id])+distribution_id;
                long double prediction_error=bin_id*fdwz_prediction_error_bin_width;
                output<<field<<','<<params.rel_eb[field]<<','<<params.abs_eb[field]<<','<<tree_id<<','<<group<<','<<bin_id<<','<<fdwz_prediction_error_bin_width<<','<<prediction_error<<','<<count<<'\n';
            }
        }
        printf("FDWZ tree prediction-error distribution written to %s\n",params.fdwz_tree_prediction_error_distribution_path.c_str());
    }

    // 执行单个field的完整FDWZ压缩；模板参数决定是否编译计时和消融统计代码。
    template<bool measure_phase_time,bool analyze_tree_distributions,bool analyze_prediction_locality,typename T>
    void fdwz_compress_impl(mrcz::param_class<T> &params,string field,uint32_t huffman_tree_num)
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
        uint32_t n012=n0*n1*n2;
        uint32_t corner_num=1;
        for (uint32_t dim=0;dim<3;dim++)
            if (header.c_element_shape[dim]>1)
                corner_num*=2;
        uint32_t residual_num=n012-corner_num;
        T abs_eb=params.abs_eb[field];
        T inv_2_eb=1/(2*abs_eb);
        if(params.fdwz_sample_interval==0)
            throw runtime_error("FDWZ sampling interval必须为正整数");
        // 与SZ3计时口径一致，总压缩时间包含内部工作区分配、plan生成、决策、预测量化和无损压缩。
        double start_time=omp_get_wtime();

        // 三个坐标方向的hierarchy和六种维度顺序只依赖element shape与reference-space GLL坐标。
        auto plan_set=fdwz_build_plan_set(n0,n1,n2,x0,x1,x2,true);
        auto &corner_positions=plan_set.corner_positions;
        auto &linear_plans=plan_set.linear_plans;
        auto &cubic_plans=plan_set.cubic_plans;
        auto all_orders=fdwz_all_dimension_orders();
        bool residual_qb_range_safe=fdwz_residual_qb_range_is_safe(params.data_range_min[field],params.data_range_max[field],abs_eb,cubic_plans);
        auto &prediction_level_offsets=cubic_plans[0].level_offsets;
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
        vector<unsigned char> selected_order_ids(blocks.size(),0);
        vector<uint64_t> sample_element_counts(blocks.size(),0);
        vector<uint64_t> unpredictable_counts(blocks.size(),0);
        vector<vector<fdwz_element_index_stencil<T>>> element_index_stencils;
        if constexpr(analyze_prediction_locality)
        {
            element_index_stencils.resize(blocks.size());
            for(uint64_t block_id=0;block_id<blocks.size();block_id++)
                element_index_stencils[block_id]=fdwz_generate_element_index_stencils<T>(blocks[block_id].element_num);
        }
        double element_local_absolute_error_sum=0,element_local_squared_error_sum=0;
        double element_index_absolute_error_sum=0,element_index_squared_error_sum=0;
        uint64_t element_local_value_num=0,element_index_value_num=0;
        vector<unique_ptr<int32_t[]>> partition_prediction_error_bins(analyze_tree_distributions?blocks.size():0);
        array<uint32_t,5> tree_value_num{};
        if (huffman_tree_num==1)
            tree_value_num[0]=n012;
        else
        {
            tree_value_num[0]=corner_num;
            for (uint32_t group_id=0;group_id+1<residual_group_offsets.size();group_id++)
                tree_value_num[group_id+1]=residual_group_offsets[group_id+1]-residual_group_offsets[group_id];
        }
        if constexpr(analyze_tree_distributions)
            for(uint64_t block_id=0;block_id<blocks.size();block_id++)
                partition_prediction_error_bins[block_id]=make_unique_for_overwrite<int32_t[]>(n012*blocks[block_id].element_num);
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
                auto &block=blocks[block_id];
                uint64_t element_num=block.end_index-block.begin_index;
                uint64_t sample_num=min<uint64_t>(element_num,max<uint64_t>(params.fdwz_min_sample_num,1+(element_num-1)/params.fdwz_sample_interval));
                sample_element_counts[block_id]=sample_num;
                auto &original_sample=sampling_workspace.original_sample;
                original_sample.resize(sample_num*n012);
                for (uint64_t sample_id=0;sample_id<sample_num;sample_id++)
                {
                    uint64_t relative_element_id=((2*sample_id+1)*element_num)/(2*sample_num);
                    memcpy(original_sample.data()+sample_id*n012,data.data()+(block.begin_index+relative_element_id)*n012,n012*sizeof(T));
                }
                uint32_t best_plan_id=0;
                uint64_t best_size=fdwz_sample_compressed_size(original_sample,n012,corner_positions,linear_plans[best_plan_id],residual_group_offsets,huffman_tree_num,abs_eb,best_plan_id,sampling_workspace);
                for(uint32_t plan_id=1;plan_id<6;plan_id++)
                {
                    uint64_t candidate_size=fdwz_sample_compressed_size(original_sample,n012,corner_positions,linear_plans[plan_id],residual_group_offsets,huffman_tree_num,abs_eb,plan_id,sampling_workspace);
                    if(candidate_size<best_size)
                        best_size=candidate_size,best_plan_id=plan_id;
                }
                selected_order_ids[block_id]=best_plan_id;
                if constexpr(measure_phase_time)
                    decision_time[block_id]=omp_get_wtime()-block_decision_start_time;
            }
            #pragma omp master
            {
                if constexpr(measure_phase_time)
                    decision_wall_time=omp_get_wtime()-decision_start_time;
            }
            #pragma omp barrier
            if constexpr(analyze_prediction_locality)
            {
                // 两种方法均使用量化前的原始值，并限到相同的目标集合：全部非pivot element中的非角点。
                #pragma omp for schedule(static) reduction(+:element_local_absolute_error_sum,element_local_squared_error_sum,element_local_value_num)
                for(uint64_t block_id=0;block_id<blocks.size();block_id++)
                {
                    auto &block=blocks[block_id];
                    auto &element_stencils=element_index_stencils[block_id];
                    auto &prediction_stencils=cubic_plans[selected_order_ids[block_id]].stencils;
                    for(const auto &element_stencil:element_stencils)
                    {
                        const T *element=data.data()+(block.begin_index+element_stencil.target)*n012;
                        for(const auto &stencil:prediction_stencils)
                        {
                            double prediction=0;
                            for(uint32_t reference_id=0;reference_id<stencil.point_num;reference_id++)
                                prediction+=double(stencil.weights[reference_id])*double(element[stencil.indices[reference_id]]);
                            double error=double(element[stencil.target])-prediction;
                            element_local_absolute_error_sum+=fabs(error);
                            element_local_squared_error_sum+=error*error;
                            element_local_value_num++;
                        }
                    }
                }
                #pragma omp for schedule(static) reduction(+:element_index_absolute_error_sum,element_index_squared_error_sum,element_index_value_num)
                for(uint64_t block_id=0;block_id<blocks.size();block_id++)
                {
                    auto &block=blocks[block_id];
                    auto &element_stencils=element_index_stencils[block_id];
                    // 所有dimension order覆盖相同的非角点集合，因此plan 0仅用于列举element-local位置。
                    auto &local_targets=cubic_plans[0].stencils;
                    for(const auto &element_stencil:element_stencils)
                    {
                        const T *target_element=data.data()+(block.begin_index+element_stencil.target)*n012;
                        for(const auto &local_target:local_targets)
                        {
                            double prediction=0;
                            for(uint32_t reference_id=0;reference_id<element_stencil.point_num;reference_id++)
                            {
                                const T *reference_element=data.data()+(block.begin_index+element_stencil.indices[reference_id])*n012;
                                prediction+=double(element_stencil.weights[reference_id])*double(reference_element[local_target.target]);
                            }
                            double error=double(target_element[local_target.target])-prediction;
                            element_index_absolute_error_sum+=fabs(error);
                            element_index_squared_error_sum+=error*error;
                            element_index_value_num++;
                        }
                    }
                }
                #pragma omp master
                {
                    double field_range=double(params.data_range_max[field])-double(params.data_range_min[field]);
                    double element_local_nrmse=sqrt(element_local_squared_error_sum/element_local_value_num)/field_range;
                    double element_local_nmae=(element_local_absolute_error_sum/element_local_value_num)/field_range;
                    double element_index_nrmse=sqrt(element_index_squared_error_sum/element_index_value_num)/field_range;
                    double element_index_nmae=(element_index_absolute_error_sum/element_index_value_num)/field_range;
                    printf("FDWZ prediction locality: element-index pivot stride=%u\n",fdwz_element_index_pivot_stride);
                    printf("FDWZ prediction locality: method=element-local, NRMSE=%.10e, NMAE=%.10e, value num=%lu\n",element_local_nrmse,element_local_nmae,element_local_value_num);
                    printf("FDWZ prediction locality: method=element-index, NRMSE=%.10e, NMAE=%.10e, value num=%lu\n",element_index_nrmse,element_index_nmae,element_index_value_num);
                }
                #pragma omp barrier
            }
            // partition覆盖互不重叠的element和QB范围，因此各线程可以安全地原地覆盖data并写入各自的bitstream。
            #pragma omp for schedule(static)
            for (uint64_t block_id=0;block_id<blocks.size();block_id++)
            {
                uint32_t selected_order_id=selected_order_ids[block_id];
                auto &prediction_stencils=cubic_plans[selected_order_id].stencils;
                unsigned char selected_plan_id=selected_order_id+6;
                int32_t *block_prediction_error_bins=partition_qb[block_id].get();
                if constexpr(analyze_tree_distributions)
                    block_prediction_error_bins=partition_prediction_error_bins[block_id].get();
                if(residual_qb_range_safe)
                    fdwz_compress_block<measure_phase_time,false,analyze_tree_distributions>(params,block_id,n012,corner_num,residual_num,corner_positions,prediction_stencils,residual_group_offsets,tree_value_num,huffman_tree_num,abs_eb,inv_2_eb,selected_plan_id,partition_qb[block_id].get(),block_prediction_error_bins,unpredictable_counts,prediction_quantization_time,huffman_time,zstd_time);
                else
                    fdwz_compress_block<measure_phase_time,true,analyze_tree_distributions>(params,block_id,n012,corner_num,residual_num,corner_positions,prediction_stencils,residual_group_offsets,tree_value_num,huffman_tree_num,abs_eb,inv_2_eb,selected_plan_id,partition_qb[block_id].get(),block_prediction_error_bins,unpredictable_counts,prediction_quantization_time,huffman_time,zstd_time);
            }
        }
        double end_time=omp_get_wtime();
        if constexpr(analyze_tree_distributions)
        {
            fdwz_write_tree_qb_distribution(params,field,huffman_tree_num,tree_value_num,partition_qb);
            fdwz_write_tree_prediction_error_distribution(params,field,huffman_tree_num,tree_value_num,partition_prediction_error_bins);
        }
        vector<uint64_t> configuration_counts(all_orders.size(),0);
        for(auto order_id:selected_order_ids)
            configuration_counts[order_id]++;
        printf("Sampled configurations:");
        for(uint32_t order_id=0;order_id<configuration_counts.size();order_id++)
            if(configuration_counts[order_id]>0)
                printf(" %c%c%c-cubic=%lu","XYZ"[all_orders[order_id][0]],"XYZ"[all_orders[order_id][1]],"XYZ"[all_orders[order_id][2]],configuration_counts[order_id]);
        printf("\n");
        printf("FDWZ QB width= 32 bits, unpredictable value num= %lu, residual QB range check= %s\n",accumulate(unpredictable_counts.begin(),unpredictable_counts.end(),uint64_t(0)),residual_qb_range_safe?"field":"per-value");
        uint64_t total_sample_element_num=accumulate(sample_element_counts.begin(),sample_element_counts.end(),uint64_t(0));
        printf("FDWZ sampling interval= %lu, minimum samples per partition= %lu, sampled element num= %lu\n",params.fdwz_sample_interval,params.fdwz_min_sample_num,total_sample_element_num);
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
        {
            huffman_tree_num=if_max_huffman_tree?0:1;
            if(if_analyze_fdwz_tree_distributions)
            {
                if(if_analyze_fdwz_prediction_locality)
                    return detail::fdwz_compress_impl<measure_phase_time,true,true>(*this,field,huffman_tree_num);
                return detail::fdwz_compress_impl<measure_phase_time,true,false>(*this,field,huffman_tree_num);
            }
            if(if_analyze_fdwz_prediction_locality)
                return detail::fdwz_compress_impl<measure_phase_time,false,true>(*this,field,huffman_tree_num);
        }
        detail::fdwz_compress_impl<measure_phase_time,false,false>(*this,field,huffman_tree_num);
    }
}
