#pragma once

#include "mrcz.hpp"
#include "fdwz_common.hpp"

namespace mrcz::detail
{
    template<bool measure_phase_time,typename T,typename QB>
    void fdwz_decompress_block(mrcz::param_class<T> &params,uint64_t block_id,uint32_t n012,uint32_t corner_num,uint32_t residual_num,vector<uint32_t> &corner_positions,
        vector<fdwz_prediction_stencil<T>> &prediction_stencils,vector<uint32_t> &residual_group_offsets,vector<uint32_t> &residual_group_ids,vector<uint32_t> &residual_group_local_ids,
        uint32_t huffman_tree_num,T abs_eb,const unsigned char *huffman_bs,const unsigned char *huffman_bs_end,vector<double> &huffman_time,vector<double> &prediction_reconstruction_time)
    {
        auto &block=params.partitions[block_id];
        auto &decompressed_data=params.decompressed_data;
        uint64_t begin_index=block.begin_index,end_index=block.end_index,element_num=end_index-begin_index;
        array<uint32_t,5> tree_value_num{};
        if(huffman_tree_num==1)
            tree_value_num[0]=n012;
        else
        {
            tree_value_num[0]=corner_num;
            for(uint32_t group_id=0;group_id+1<residual_group_offsets.size();group_id++)
                tree_value_num[group_id+1]=residual_group_offsets[group_id+1]-residual_group_offsets[group_id];
        }
        array<uint64_t,6> tree_offsets{};
        for(uint32_t tree_id=0;tree_id<huffman_tree_num;tree_id++)
            tree_offsets[tree_id+1]=tree_offsets[tree_id]+tree_value_num[tree_id]*element_num;
        double phase_start_time=0;
        if constexpr(measure_phase_time)
            phase_start_time=omp_get_wtime();
        auto block_qb=make_unique_for_overwrite<QB[]>(n012*element_num);
        for(uint32_t tree_id=0;tree_id<huffman_tree_num;tree_id++)
        {
            mrcz::huffman_encoder_class<QB> encoder;
            encoder.load(huffman_bs);
            encoder.decode(huffman_bs,block_qb.get()+tree_offsets[tree_id],tree_offsets[tree_id+1]-tree_offsets[tree_id]);
        }
        if constexpr(measure_phase_time)
            huffman_time[block_id]=omp_get_wtime()-phase_start_time;
        const unsigned char *unpredictable_values=huffman_bs;
        auto corner_qb_at=[&](uint64_t relative_element_id,uint32_t corner_id)->QB
        {
            return huffman_tree_num==1?block_qb[relative_element_id*n012+corner_id]:block_qb[relative_element_id*corner_num+corner_id];
        };
        auto residual_qb_at=[&](uint64_t relative_element_id,uint32_t residual_id)->QB
        {
            if(huffman_tree_num==1)
                return block_qb[relative_element_id*n012+corner_num+residual_id];
            uint32_t group_id=residual_group_ids[residual_id];
            uint32_t group_size=residual_group_offsets[group_id+1]-residual_group_offsets[group_id];
            return block_qb[tree_offsets[group_id+1]+relative_element_id*group_size+residual_group_local_ids[residual_id]];
        };
        if constexpr(measure_phase_time)
            phase_start_time=omp_get_wtime();
        auto recover_value=[&](T prediction,QB qb)->T
        {
            if constexpr(is_same_v<QB,int32_t>)
                if(qb==numeric_limits<int32_t>::min())
                {
                    if(unpredictable_values+sizeof(T)>huffman_bs_end)
                        throw runtime_error("FDWZ unpredictable原值数量不足");
                    T value;
                    memcpy(&value,unpredictable_values,sizeof(T));
                    unpredictable_values+=sizeof(T);
                    return value;
                }
            return prediction-2*abs_eb*qb;
        };
        for(uint64_t i=begin_index;i<end_index;i++)
        {
            uint64_t relative_element_id=i-begin_index;
            for(uint32_t corner_id=0;corner_id<corner_num;corner_id++)
            {
                uint32_t target=corner_positions[corner_id];
                decompressed_data[i*n012+target]=recover_value(0,corner_qb_at(relative_element_id,corner_id));
            }
        }
        if(huffman_tree_num==1)
        {
            for(uint64_t i=begin_index;i<end_index;i++)
            {
                T *element=decompressed_data.data()+i*n012;
                uint64_t relative_element_id=i-begin_index;
                for(uint32_t residual_id=0;residual_id<residual_num;residual_id++)
                {
                    auto &stencil=prediction_stencils[residual_id];
                    T prediction=fdwz_predict_from_stencil(element,stencil);
                    element[stencil.target]=recover_value(prediction,residual_qb_at(relative_element_id,residual_id));
                }
            }
        }
        else
        {
            constexpr uint64_t batch_size=4;
            uint64_t batch_end=begin_index+(element_num/batch_size)*batch_size;
            for(uint64_t batch_begin=begin_index;batch_begin<batch_end;batch_begin+=batch_size)
                for(uint32_t group_id=0;group_id+1<residual_group_offsets.size();group_id++)
                {
                    uint32_t group_begin=residual_group_offsets[group_id],group_size=residual_group_offsets[group_id+1]-group_begin;
                    for(uint32_t local_id=0;local_id<group_size;local_id++)
                    {
                        uint32_t residual_id=group_begin+local_id;
                        auto &stencil=prediction_stencils[residual_id];
                        #pragma GCC unroll 4
                        for(uint64_t lane=0;lane<batch_size;lane++)
                        {
                            uint64_t i=batch_begin+lane,relative_element_id=i-begin_index;
                            T *element=decompressed_data.data()+i*n012;
                            T prediction=fdwz_predict_from_stencil(element,stencil);
                            element[stencil.target]=recover_value(prediction,residual_qb_at(relative_element_id,residual_id));
                        }
                    }
                }
            for(uint64_t i=batch_end;i<end_index;i++)
            {
                T *element=decompressed_data.data()+i*n012;
                uint64_t relative_element_id=i-begin_index;
                for(uint32_t group_id=0;group_id+1<residual_group_offsets.size();group_id++)
                {
                    uint32_t group_begin=residual_group_offsets[group_id],group_size=residual_group_offsets[group_id+1]-group_begin;
                    for(uint32_t local_id=0;local_id<group_size;local_id++)
                    {
                        uint32_t residual_id=group_begin+local_id;
                        auto &stencil=prediction_stencils[residual_id];
                        T prediction=fdwz_predict_from_stencil(element,stencil);
                        element[stencil.target]=recover_value(prediction,residual_qb_at(relative_element_id,residual_id));
                    }
                }
            }
        }
        if(unpredictable_values!=huffman_bs_end)
            throw runtime_error("FDWZ unpredictable原值数量过多");
        if constexpr(measure_phase_time)
            prediction_reconstruction_time[block_id]=omp_get_wtime()-phase_start_time;
    }

    // 执行单个field的完整FDWZ解压，并使用与压缩端相同的shape和配置重建全部预测plan。
    template<bool measure_phase_time,typename T>
    void fdwz_decompress_impl(mrcz::param_class<T> &params,string field,uint32_t huffman_tree_num,bool use_gll_weights,bool use_gll_order)
    {
        auto &header=params.header;
        auto &blocks=params.partitions;
        auto &decompressed_data=params.decompressed_data;
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
        decompressed_data.resize(params.element_num_to_process*n012);
        // 与其他backend一致，输出数组分配不计时；总解压时间包含hierarchy、plan生成、无损解压和预测重构。
        double start_time=omp_get_wtime();

        // 解压端只重建正式预测使用的六种cubic plans，linear proxy只在压缩端选择维度顺序时使用。
        auto plan_set=fdwz_build_plan_set(n0,n1,n2,x0,x1,x2,false,use_gll_weights,use_gll_order);
        auto &corner_positions=plan_set.corner_positions;
        auto &cubic_plans=plan_set.cubic_plans;
        auto &prediction_level_offsets=cubic_plans[0].level_offsets;
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
        vector<uint32_t> residual_group_ids(residual_num),residual_group_local_ids(residual_num);
        for(uint32_t group_id=0;group_id+1<residual_group_offsets.size();group_id++)
            for (uint32_t residual_id=residual_group_offsets[group_id];residual_id<residual_group_offsets[group_id+1];residual_id++)
                residual_group_ids[residual_id]=group_id,residual_group_local_ids[residual_id]=residual_id-residual_group_offsets[group_id];

        vector<double> zstd_time(measure_phase_time?blocks.size():0,0),huffman_time(measure_phase_time?blocks.size():0,0),prediction_reconstruction_time(measure_phase_time?blocks.size():0,0);
        // 每个partition拥有独立bitstream和不重叠的输出范围，可以并行完成解码与重构。
        #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
        for (uint64_t block_id=0;block_id<blocks.size();block_id++)
        {
            auto &block=blocks[block_id];
            // 先用Zstd恢复Huffman payload，再读取首字节中的plan id。
            double phase_start_time=0;
            if constexpr(measure_phase_time)
                phase_start_time=omp_get_wtime();
            uint64_t huffman_bs_len=ZSTD_getFrameContentSize(block.zstd_bs.data(),block.zstd_bs.size());
            vector<unsigned char> huffman_bs_vec(huffman_bs_len);
            if (ZSTD_decompress(huffman_bs_vec.data(),huffman_bs_len,block.zstd_bs.data(),block.zstd_bs.size())!=huffman_bs_len)
                throw runtime_error("ZSTD_decompress error in decompress()");
            if constexpr(measure_phase_time)
                zstd_time[block_id]=omp_get_wtime()-phase_start_time;
            const unsigned char *huffman_bs=huffman_bs_vec.data();
            if (huffman_bs_len<1)
                throw runtime_error("FDWZ采样配置的bitstream过短");
            uint32_t selected_plan_id=*huffman_bs++;
            if(selected_plan_id<6 or selected_plan_id>=6+cubic_plans.size())
                throw runtime_error("FDWZ采样配置metadata无效");
            auto &prediction_stencils=cubic_plans[selected_plan_id-6].stencils;
            fdwz_decompress_block<measure_phase_time,T,int32_t>(params,block_id,n012,corner_num,residual_num,corner_positions,prediction_stencils,residual_group_offsets,residual_group_ids,residual_group_local_ids,huffman_tree_num,abs_eb,huffman_bs,huffman_bs_vec.data()+huffman_bs_len,huffman_time,prediction_reconstruction_time);
        }
        double end_time=omp_get_wtime();
        if constexpr(measure_phase_time)
        {
            double total_zstd_time=0,total_huffman_time=0,total_prediction_reconstruction_time=0;
            for(uint64_t block_id=0;block_id<blocks.size();block_id++)
            {
                total_zstd_time+=zstd_time[block_id];
                total_huffman_time+=huffman_time[block_id];
                total_prediction_reconstruction_time+=prediction_reconstruction_time[block_id];
            }
            double measured_stage_time=total_zstd_time+total_huffman_time+total_prediction_reconstruction_time;
            printf("FDWZ decompression accumulated partition CPU times: Zstd= %f, Huffman= %f, QB consumption+prediction+reconstruction= %f seconds\n",total_zstd_time,total_huffman_time,total_prediction_reconstruction_time);
            if(measured_stage_time>0)
                printf("FDWZ decompression stage shares: Zstd= %.2f%%, Huffman= %.2f%%, QB consumption+prediction+reconstruction= %.2f%%\n",100*total_zstd_time/measured_stage_time,100*total_huffman_time/measured_stage_time,100*total_prediction_reconstruction_time/measured_stage_time);
        }
        printf("Decompression time= %f seconds\n",end_time-start_time);
        printf("\033[91mField= %s, decompression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),params.element_num_to_process*n012*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
    }
}

namespace mrcz
{
    // 根据正式版或消融版配置选择Huffman tree数量，并进入FDWZ解压实现。
    template<typename T>
    template<bool measure_phase_time,bool ablation>
    void param_class<T>::decompress_fdwz(string field)
    {
        uint32_t huffman_tree_num=0;
        if constexpr(ablation)
        {
            huffman_tree_num=fdwz_huffman_tree_num>0?fdwz_huffman_tree_num:(if_max_huffman_tree?0:1);
            return detail::fdwz_decompress_impl<measure_phase_time>(*this,field,huffman_tree_num,fdwz_use_gll_weights,fdwz_use_gll_order);
        }
        detail::fdwz_decompress_impl<measure_phase_time>(*this,field,huffman_tree_num,true,true);
    }
}
