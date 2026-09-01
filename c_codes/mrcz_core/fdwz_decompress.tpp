#pragma once


#include "mrcz.hpp"
#include "fdwz_compress.tpp"

namespace mrcz::detail
{
    // 执行单个field的完整FDWZ解压，并使用与压缩端相同的shape和配置重建全部预测plan。
    template<bool measure_phase_time,typename T>
    void fdwz_decompress_impl(mrcz::param_class<T> &params,string field,uint32_t huffman_tree_num)
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
        uint32_t n12=n1*n2;
        uint32_t n012=n0*n1*n2;
        uint32_t corner_num=1;
        for (uint32_t dim=0;dim<3;dim++)
            if (header.c_element_shape[dim]>1)
                corner_num*=2;
        uint32_t residual_num=n012-corner_num;
        T abs_eb=params.abs_eb[field];
        // 总解压时间包含hierarchy、plan生成、无损解压和预测重构，与压缩端使用相同的计时边界。
        double start_time=omp_get_wtime();
        decompressed_data.resize(params.element_num_to_process*n012);

        // 解压端根据坐标和shape确定性地重建全部候选plan，bitstream只需保存每个partition选择的plan id。
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
        auto &prediction_level_offsets=sampled_plans[0].level_offsets;
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
        int64_t *corner_qb=huffman_tree_num==1?nullptr:new int64_t[corner_num*params.element_num_to_process];
        vector<int64_t*> residual_group_qb;
        for (uint32_t group_id=0;group_id+1<residual_group_offsets.size();group_id++)
            residual_group_qb.push_back(new int64_t[params.element_num_to_process*(residual_group_offsets[group_id+1]-residual_group_offsets[group_id])]);
        vector<uint32_t> residual_group_ids(residual_num),residual_group_local_ids(residual_num);
        for (uint32_t group_id=0;group_id<residual_group_qb.size();group_id++)
            for (uint32_t residual_id=residual_group_offsets[group_id];residual_id<residual_group_offsets[group_id+1];residual_id++)
                residual_group_ids[residual_id]=group_id,residual_group_local_ids[residual_id]=residual_id-residual_group_offsets[group_id];

        vector<double> zstd_time(measure_phase_time?blocks.size():0,0),huffman_time(measure_phase_time?blocks.size():0,0),prediction_reconstruction_time(measure_phase_time?blocks.size():0,0);
        // 每个partition拥有独立bitstream和不重叠的输出范围，可以并行完成解码与重构。
        #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
        for (uint64_t block_id=0;block_id<blocks.size();block_id++)
        {
            auto &block=blocks[block_id];
            uint64_t begin_index=block.begin_index;
            uint64_t end_index=block.end_index;
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
            if (selected_plan_id>=sampled_plans.size())
                throw runtime_error("FDWZ采样配置metadata无效");
            auto &block_prediction_stencils=sampled_plans[selected_plan_id].stencils;
            uint64_t element_num=end_index-begin_index;
            int64_t *block_corner_qb=huffman_tree_num==1?nullptr:corner_qb+corner_num*begin_index;
            vector<int64_t*> block_residual_group_qb;
            for (uint32_t group_id=0;group_id<residual_group_qb.size();group_id++)
                block_residual_group_qb.push_back(residual_group_qb[group_id]+begin_index*(residual_group_offsets[group_id+1]-residual_group_offsets[group_id]));
            vector<int64_t> block_qb;
            if constexpr(measure_phase_time)
                phase_start_time=omp_get_wtime();
            if (huffman_tree_num==1)
            {
                block_qb.resize(n012*element_num);
                mrcz::huffman_encoder_class<int64_t> encoder;
                encoder.load(huffman_bs);
                encoder.decode(huffman_bs,block_qb.data(),block_qb.size());
            }
            else
            {
                mrcz::huffman_encoder_class<int64_t> corner_encoder;
                corner_encoder.load(huffman_bs);
                corner_encoder.decode(huffman_bs,block_corner_qb,corner_num*element_num);
                for (uint32_t group_id=0;group_id<block_residual_group_qb.size();group_id++)
                {
                    uint64_t group_value_num=(residual_group_offsets[group_id+1]-residual_group_offsets[group_id])*element_num;
                    mrcz::huffman_encoder_class<int64_t> encoder;
                    encoder.load(huffman_bs);
                    encoder.decode(huffman_bs,block_residual_group_qb[group_id],group_value_num);
                }
            }
            if constexpr(measure_phase_time)
                huffman_time[block_id]=omp_get_wtime()-phase_start_time;
            // 两个访问器隐藏一棵tree的element连续布局和多棵tree的分组布局差异。
            auto corner_qb_at=[&](uint64_t element_id,uint32_t corner_id)->int64_t
            {
                if (huffman_tree_num==1)
                    return block_qb[(element_id-begin_index)*n012+corner_id];
                return corner_qb[element_id*corner_num+corner_id];
            };
            auto residual_qb_at=[&](uint64_t element_id,uint32_t residual_id)->int64_t
            {
                if (huffman_tree_num==1)
                    return block_qb[(element_id-begin_index)*n012+corner_num+residual_id];
                uint32_t group_id=residual_group_ids[residual_id];
                uint32_t group_size=residual_group_offsets[group_id+1]-residual_group_offsets[group_id];
                return residual_group_qb[group_id][element_id*group_size+residual_group_local_ids[residual_id]];
            };
            if constexpr(measure_phase_time)
                phase_start_time=omp_get_wtime();
            // 先恢复角点，再按plan顺序恢复非角点，保证每个stencil引用的值都已经重构。
            for (uint64_t i=begin_index;i<end_index;i++)
                for (uint32_t corner_id=0;corner_id<corner_num;corner_id++)
                    decompressed_data[i*n012+corner_positions[corner_id]]=corner_qb_at(i,corner_id)*2*abs_eb;
            for (uint64_t i=begin_index;i<end_index;i++)
            {
                T *element=decompressed_data.data()+i*n012;
                uint32_t residual_id=0;
                for (auto &stencil:block_prediction_stencils)
                {
                    T prediction=fdwz_predict_from_stencil(element,stencil);
                    element[stencil.target]=prediction-2*abs_eb*residual_qb_at(i,residual_id++);
                }
                if (residual_id!=residual_num)
                    throw runtime_error("FDWZ重构读取的QB数量错误");
            }
            if constexpr(measure_phase_time)
                prediction_reconstruction_time[block_id]=omp_get_wtime()-phase_start_time;
        }
        double end_time=omp_get_wtime();
        delete[] corner_qb;
        for (auto group_qb:residual_group_qb)
            delete[] group_qb;
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
            huffman_tree_num=if_max_huffman_tree?0:1;
        detail::fdwz_decompress_impl<measure_phase_time>(*this,field,huffman_tree_num);
    }
}
