#pragma once


#include "mrcz.hpp"
#include "fdwz_compress.tpp"

template<typename T>
void sz3_3d_gll_decompress_impl(mrcz::param_class<T> &params,string field,uint32_t huffman_tree_num,vector<uint32_t> dimension_order={0,1,2},bool use_cubic=false,bool use_surface_first=false,bool use_normal_boundary=false,bool use_index_order=false,bool use_index_midpoint=false,T reference_eb_scale=1,bool use_position_bias=false,bool sample_configuration=false)
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
    uint32_t seed_num=use_index_order?1:corner_num;
    uint32_t residual_num=n012-seed_num;
    T abs_eb=params.abs_eb[field];
    decompressed_data.resize(params.element_num_to_process*n012);
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
    if (use_position_bias and corner_positions.size()!=seed_num)
        throw runtime_error("SZ3位置bias的角点数量错误");
    T reference_abs_eb=abs_eb*reference_eb_scale;
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
    int64_t *corner_qb=huffman_tree_num==1?nullptr:new int64_t[seed_num*params.element_num_to_process];
    int64_t *no_corner_err_qb=huffman_tree_num==2?new int64_t[params.element_num_to_process*residual_num]:nullptr;
    vector<int64_t*> residual_group_qb;
    for (uint32_t group_id=0;group_id+1<residual_group_offsets.size();group_id++)
        residual_group_qb.push_back(new int64_t[params.element_num_to_process*(residual_group_offsets[group_id+1]-residual_group_offsets[group_id])]);
    vector<uint32_t> residual_group_ids(residual_num),residual_group_local_ids(residual_num);
    for (uint32_t group_id=0;group_id<residual_group_qb.size();group_id++)
        for (uint32_t residual_id=residual_group_offsets[group_id];residual_id<residual_group_offsets[group_id+1];residual_id++)
            residual_group_ids[residual_id]=group_id,residual_group_local_ids[residual_id]=residual_id-residual_group_offsets[group_id];

    double start_time=omp_get_wtime();
    double zstd_time=0,huffman_time=0,prediction_reconstruction_time=0;
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus) reduction(+:zstd_time,huffman_time,prediction_reconstruction_time)
    for (uint64_t block_id=0;block_id<blocks.size();block_id++)
    {
        auto &block=blocks[block_id];
        uint64_t begin_index=block.begin_index;
        uint64_t end_index=block.end_index;
        // 先无损解压
        double zstd_start_time=omp_get_wtime();
        uint64_t huffman_bs_len=ZSTD_getFrameContentSize(block.zstd_bs.data(),block.zstd_bs.size());
        vector<unsigned char> huffman_bs_vec(huffman_bs_len);
        if (ZSTD_decompress(huffman_bs_vec.data(),huffman_bs_len,block.zstd_bs.data(),block.zstd_bs.size())!=huffman_bs_len)
            throw runtime_error("ZSTD_decompress error in decompress()");
        zstd_time+=omp_get_wtime()-zstd_start_time;
        const unsigned char *huffman_bs=huffman_bs_vec.data();
        uint32_t selected_plan_id=0;
        if (sample_configuration)
        {
            if (huffman_bs_len<1)
                throw runtime_error("SZ3采样配置的bitstream过短");
            selected_plan_id=*huffman_bs++;
            if (selected_plan_id>=sampled_plans.size())
                throw runtime_error("SZ3采样配置metadata无效");
        }
        vector<sz3_prediction_stencil<T>> *block_prediction_stencils=sample_configuration?&sampled_plans[selected_plan_id].stencils:nullptr;
        vector<int32_t> position_bias(n012,0);
        T bias_unit=abs_eb/8;
        if (use_position_bias)
        {
            uint64_t position_bias_bytes=n012*sizeof(int32_t);
            if (huffman_bs_len<position_bias_bytes)
                throw runtime_error("SZ3位置bias的bitstream过短");
            memcpy(position_bias.data(),huffman_bs,position_bias_bytes);
            huffman_bs+=position_bias_bytes;
        }
        uint64_t element_num=end_index-begin_index;
        int64_t *block_corner_qb=huffman_tree_num==1?nullptr:corner_qb+seed_num*begin_index;
        int64_t *block_no_corner_err_qb=huffman_tree_num==2?no_corner_err_qb+residual_num*begin_index:nullptr;
        vector<int64_t*> block_residual_group_qb;
        for (uint32_t group_id=0;group_id<residual_group_qb.size();group_id++)
            block_residual_group_qb.push_back(residual_group_qb[group_id]+begin_index*(residual_group_offsets[group_id+1]-residual_group_offsets[group_id]));
        vector<int64_t> block_qb;
        double huffman_start_time=omp_get_wtime();
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
            corner_encoder.decode(huffman_bs,block_corner_qb,seed_num*element_num);
            if (huffman_tree_num==2)
            {
                mrcz::huffman_encoder_class<int64_t> err_encoder;
                err_encoder.load(huffman_bs);
                err_encoder.decode(huffman_bs,block_no_corner_err_qb,residual_num*element_num);
            }
            else
            {
                for (uint32_t group_id=0;group_id<block_residual_group_qb.size();group_id++)
                {
                    uint64_t group_value_num=(residual_group_offsets[group_id+1]-residual_group_offsets[group_id])*element_num;
                    mrcz::huffman_encoder_class<int64_t> encoder;
                    encoder.load(huffman_bs);
                    encoder.decode(huffman_bs,block_residual_group_qb[group_id],group_value_num);
                }
            }
        }
        huffman_time+=omp_get_wtime()-huffman_start_time;
        auto corner_qb_at=[&](uint64_t element_id,uint32_t corner_id)->int64_t
        {
            if (huffman_tree_num==1)
                return block_qb[(element_id-begin_index)*n012+corner_id];
            return corner_qb[element_id*seed_num+corner_id];
        };
        auto residual_qb_at=[&](uint64_t element_id,uint32_t residual_id)->int64_t
        {
            if (huffman_tree_num==1)
                return block_qb[(element_id-begin_index)*n012+seed_num+residual_id];
            if (huffman_tree_num==2)
                return no_corner_err_qb[element_id*residual_num+residual_id];
            uint32_t group_id=residual_group_ids[residual_id];
            uint32_t group_size=residual_group_offsets[group_id+1]-residual_group_offsets[group_id];
            return residual_group_qb[group_id][element_id*group_size+residual_group_local_ids[residual_id]];
        };
        double prediction_reconstruction_start_time=omp_get_wtime();
        // 重构每个element角上的数据
        if (use_position_bias)
            for (uint64_t i=begin_index;i<end_index;i++)
                for (uint32_t corner_id=0;corner_id<seed_num;corner_id++)
                {
                    uint32_t target=corner_positions[corner_id];
                    decompressed_data[i*n012+target]=position_bias[target]*bias_unit-2*abs_eb*corner_qb_at(i,corner_id);
                }
        else if (reference_eb_scale<1)
            for (uint64_t i=begin_index;i<end_index;i++)
                for (uint32_t corner_id=0;corner_id<seed_num;corner_id++)
                    decompressed_data[i*n012+corner_positions[corner_id]]=corner_qb_at(i,corner_id)*2*reference_abs_eb;
        else
            for (uint64_t i=begin_index;i<end_index;i++)
                for (uint32_t corner_id=0;corner_id<seed_num;corner_id++)
                    decompressed_data[i*n012+corner_positions[corner_id]]=corner_qb_at(i,corner_id)*2*abs_eb;
        for (uint64_t i=begin_index;i<end_index;i++)
        {
            T *element=decompressed_data.data()+i*n012;
            uint32_t residual_id=0;
            auto reconstruct_prediction=[&](T pred_value,uint32_t target)
            {
                if (use_position_bias)
                    pred_value+=position_bias[target]*bias_unit;
                T point_abs_eb=reference_eb_scale<1 and reference_mask[target]?reference_abs_eb:abs_eb;
                element[target]=pred_value-2*point_abs_eb*residual_qb_at(i,residual_id++);
            };
            if (sample_configuration)
            {
                for (auto &stencil:*block_prediction_stencils)
                    reconstruct_prediction(sz3_predict_from_stencil(element,stencil),stencil.target);
            }
            else if (use_index_order)
            {
                for (auto &stencil:index_order_stencils)
                {
                    T pred_value=sz3_index_order_predict(element,stencil);
                    reconstruct_prediction(pred_value,stencil.target);
                }
            }
            else if (use_normal_boundary)
            {
                for (auto &stencil:normal_boundary_stencils)
                {
                    T pred_value=sz3_normal_boundary_predict(element,stencil);
                    reconstruct_prediction(pred_value,stencil.target);
                }
            }
            else if (use_surface_first)
            {
                for (auto &stencil:surface_first_stencils)
                {
                    T pred_value=sz3_surface_first_predict(element,stencil);
                    reconstruct_prediction(pred_value,stencil.target);
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
                                        reconstruct_prediction(pred_value,mid_pos);
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
                                        reconstruct_prediction(pred_value,mid_pos);
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
                                        reconstruct_prediction(pred_value,mid_pos);
                                    }
                            level2++;
                        }
                    }
                }
            }
            if (residual_id!=residual_num)
                throw runtime_error("SZ3重构读取的QB数量错误");
        }
        prediction_reconstruction_time+=omp_get_wtime()-prediction_reconstruction_start_time;
    }
    double end_time=omp_get_wtime();
    delete[] corner_qb;
    delete[] no_corner_err_qb;
    for (auto group_qb:residual_group_qb)
        delete[] group_qb;
    double measured_stage_time=zstd_time+huffman_time+prediction_reconstruction_time;
    printf("Decompression accumulated partition times: Zstd= %f, Huffman= %f, QB consumption+prediction+reconstruction= %f seconds\n",zstd_time,huffman_time,prediction_reconstruction_time);
    if (measured_stage_time>0)
        printf("Decompression stage shares: Zstd= %.2f%%, Huffman= %.2f%%, QB consumption+prediction+reconstruction= %.2f%%\n",100*zstd_time/measured_stage_time,100*huffman_time/measured_stage_time,100*prediction_reconstruction_time/measured_stage_time);
    printf("Decompression time= %f seconds\n",end_time-start_time);
    printf("\033[91mField= %s, decompression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),params.element_num_to_process*n012*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
}

template<typename T>
void sz3_3d_gll_decompress(mrcz::param_class<T> &params,string field)
{
    sz3_3d_gll_decompress_impl(params,field,1);
}

template<typename T>
void fdwz_decompress(mrcz::param_class<T> &params,string field)
{
    sz3_3d_gll_decompress_impl(params,field,0,{0,1,2},false,false,false,false,false,T(1),false,true);
}

template<typename T>
void sz3_2_decompress(mrcz::param_class<T> &params,string field)
{
    sz3_3d_gll_decompress_impl(params,field,2);
}

template<typename T>
void sz3_3d_gll_experimental_decompress(mrcz::param_class<T> &params,string field)
{
    sz3_3d_gll_decompress_impl(params,field,5);
}

template<typename T>
void sz3_3d_gll_sampled_configuration_decompress(mrcz::param_class<T> &params,string field)
{
    auto configuration=sz3_sampled_configuration_from_method(params.method);
    sz3_3d_gll_decompress_impl(params,field,configuration.huffman_tree_num,configuration.dimension_order,configuration.use_cubic);
}

template<typename T>
void sz3_3d_gll_position_bias_decompress(mrcz::param_class<T> &params,string field)
{
    auto configuration=sz3_position_bias_configuration_from_method(params.method);
    sz3_3d_gll_decompress_impl(params,field,configuration.huffman_tree_num,configuration.dimension_order,configuration.use_cubic,false,false,false,false,T(1),true);
}

template<typename T>
void sz3_3d_gll_order_decompress(mrcz::param_class<T> &params,string field)
{
    sz3_3d_gll_decompress_impl(params,field,2,sz3_dimension_order_from_method(params.method));
}

template<typename T>
void sz3_3d_gll_cubic_decompress(mrcz::param_class<T> &params,string field)
{
    sz3_3d_gll_decompress_impl(params,field,2,{0,1,2},true);
}

template<typename T>
void sz3_3d_gll_surface_first_decompress(mrcz::param_class<T> &params,string field)
{
    sz3_3d_gll_decompress_impl(params,field,2,{0,1,2},false,true);
}

template<typename T>
void sz3_3d_gll_normal_boundary_decompress(mrcz::param_class<T> &params,string field)
{
    sz3_3d_gll_decompress_impl(params,field,2,{0,1,2},false,false,true);
}

template<typename T>
void sz3_3d_gll_index_order_decompress(mrcz::param_class<T> &params,string field)
{
    sz3_3d_gll_decompress_impl(params,field,2,{0,1,2},false,false,false,true);
}

template<typename T>
void sz3_3d_gll_index_midpoint_decompress(mrcz::param_class<T> &params,string field)
{
    sz3_3d_gll_decompress_impl(params,field,2,{0,1,2},false,false,false,false,true);
}

template<typename T>
void sz3_3d_gll_reference_eb_decompress(mrcz::param_class<T> &params,string field)
{
    sz3_3d_gll_decompress_impl(params,field,2,{0,1,2},false,false,false,false,false,T(sz3_reference_eb_scale_from_method(params.method)));
}

namespace mrcz
{
    template<typename T>
    inline void param_class<T>::reconstruct_corner_number(int64_t *corner_qb,uint64_t begin_index,uint64_t end_index,uint32_t n0,uint32_t n1,uint32_t n2,uint32_t n12,uint32_t n012,T abs_eb)
    {
        if ((n0>1) and (n1>1) and (n2>1))
            for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
            {
                T *element=decompressed_data.data()+element_id*n012;
                int64_t *element_corner_qb=corner_qb+8*element_id;
                element[0]=element_corner_qb[0]*2*abs_eb;
                element[n2-1]=element_corner_qb[1]*2*abs_eb;
                element[n12-n2]=element_corner_qb[2]*2*abs_eb;
                element[n12-1]=element_corner_qb[3]*2*abs_eb;
                element[n012-n12]=element_corner_qb[4]*2*abs_eb;
                element[n012-n12+n2-1]=element_corner_qb[5]*2*abs_eb;
                element[n012-n2]=element_corner_qb[6]*2*abs_eb;
                element[n012-1]=element_corner_qb[7]*2*abs_eb;
            }
        else if ((n0>1) and (n1>1) and (n2==1))
            for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
            {
                T *element=decompressed_data.data()+element_id*n012;
                int64_t *element_corner_qb=corner_qb+4*element_id;
                element[0]=element_corner_qb[0]*2*abs_eb;
                element[n12-n2]=element_corner_qb[1]*2*abs_eb;
                element[n012-n12]=element_corner_qb[2]*2*abs_eb;
                element[n012-n2]=element_corner_qb[3]*2*abs_eb;
            }
        else if ((n0>1) and (n1==1) and (n2>1))
            for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
            {
                T *element=decompressed_data.data()+element_id*n012;
                int64_t *element_corner_qb=corner_qb+4*element_id;
                element[0]=element_corner_qb[0]*2*abs_eb;
                element[n2-1]=element_corner_qb[1]*2*abs_eb;
                element[n012-n12]=element_corner_qb[2]*2*abs_eb;
                element[n012-n12+n2-1]=element_corner_qb[3]*2*abs_eb;
            }
        else if ((n0==1) and (n1>1) and (n2>1))
            for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
            {
                T *element=decompressed_data.data()+element_id*n012;
                int64_t *element_corner_qb=corner_qb+4*element_id;
                element[0]=element_corner_qb[0]*2*abs_eb;
                element[n2-1]=element_corner_qb[1]*2*abs_eb;
                element[n12-n2]=element_corner_qb[2]*2*abs_eb;
                element[n12-1]=element_corner_qb[3]*2*abs_eb;
            }
        else
            throw runtime_error("Error: sz3_2_decompress only supports 3D or 2D elements.\n");
    }

    template<typename T>
    inline void param_class<T>::convert_no_corner_err_eb_to_err_eb(int64_t *element_no_corner_err_qb,int64_t *element_err_qb,uint32_t n0,uint32_t n1,uint32_t n2,uint32_t n12,uint32_t n012)
    {
        if ((n0>1) and (n1>1) and (n2>1))
        {
            memcpy(element_err_qb+1,element_no_corner_err_qb,sizeof(int64_t)*(n2-2));
            memcpy(element_err_qb+n2,element_no_corner_err_qb+n2-2,sizeof(int64_t)*(n12-2*n2));
            memcpy(element_err_qb+n12-n2+1,element_no_corner_err_qb+n12-n2-2,sizeof(int64_t)*(n2-2));
            memcpy(element_err_qb+n12,element_no_corner_err_qb+n12-4,sizeof(int64_t)*(n012-2*n12));
            memcpy(element_err_qb+n012-n12+1,element_no_corner_err_qb+n012-n12-4,sizeof(int64_t)*(n2-2));
            memcpy(element_err_qb+n012-n12+n2,element_no_corner_err_qb+n012-n12+n2-6,sizeof(int64_t)*(n12-2*n2));
            memcpy(element_err_qb+n012-n2+1,element_no_corner_err_qb+n012-n2-6,sizeof(int64_t)*(n2-2));
        }
        else if ((n0>1) and (n1>1) and (n2==1))
        {
            memcpy(element_err_qb+1,element_no_corner_err_qb,sizeof(int64_t)*(n12-2));
            memcpy(element_err_qb+n12,element_no_corner_err_qb+n12-2,sizeof(int64_t)*(n012-2*n12));
            memcpy(element_err_qb+n012-n12+1,element_no_corner_err_qb+n012-n12-2,sizeof(int64_t)*(n12-2));
        }
        else if ((n0>1) and (n1==1) and (n2>1))
        {
            memcpy(element_err_qb+1,element_no_corner_err_qb,sizeof(int64_t)*(n2-2));
            memcpy(element_err_qb+n2,element_no_corner_err_qb+n2-2,sizeof(int64_t)*(n012-2*n2));
            memcpy(element_err_qb+n012-n2+1,element_no_corner_err_qb+n012-n2-2,sizeof(int64_t)*(n2-2));
        }
        else if ((n0==1) and (n1>1) and (n2>1))
        {
            memcpy(element_err_qb+1,element_no_corner_err_qb,sizeof(int64_t)*(n2-2));
            memcpy(element_err_qb+n2,element_no_corner_err_qb+n2-2,sizeof(int64_t)*(n12-2*n2));
            memcpy(element_err_qb+n12-n2+1,element_no_corner_err_qb+n12-n2-2,sizeof(int64_t)*(n2-2));
        }
        else
            throw runtime_error("Error: sz3_2_decompress only supports 3D or 2D elements.\n");
    }
}
