// XXX 这是SZ3-3D-GLL-experimental旧版12种角点流诊断实现的完整快照。
// XXX 该版本只打印Huffman+Zstd统计，不生成可解压bitstream；活动实现已改为完整I0+T格式。

#pragma once


#include <numeric>
#include <set>
#include <unordered_map>
#include "mrcz.hpp"
#include "huffman.hpp"
#include "quantize.tpp"
#include "../third_party/zstd/lib/zstd.h"

template<typename T>
void sz3_generate_weight(vector<vector<uint32_t>> &left_indices,vector<vector<uint32_t>> &right_indices,vector<vector<uint32_t>> &mid_indices,
    vector<vector<uint32_t>> &compressed_indices,vector<vector<T>> &left_weight,vector<vector<T>> &right_weight,vector<T> &x)
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
            T ideal_mid_pos=(x[left_index]+x[right_index])/2;
            for (uint32_t j=left_index+2;j<right_index;j++)
                if (fabs(x[j]-ideal_mid_pos)<fabs(x[mid_index]-ideal_mid_pos))
                    mid_index=j;
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

inline pair<double,double> sz3_quantization_statistics(int64_t *qb,uint64_t count)
{
    uint64_t zero_count=0;
    unordered_map<int64_t,uint64_t> frequencies;
    for (uint64_t i=0;i<count;i++)
    {
        zero_count+=qb[i]==0;
        frequencies[qb[i]]++;
    }
    double entropy=0;
    for (auto &item:frequencies)
    {
        double probability=static_cast<double>(item.second)/count;
        entropy-=probability*log2(probability);
    }
    return {static_cast<double>(zero_count)/count,entropy};
}

inline vector<uint32_t> sz3_corner_positions(uint32_t n0,uint32_t n1,uint32_t n2)
{
    vector<uint32_t> positions;
    vector<uint32_t> indices0={0},indices1={0},indices2={0};
    if (n0>1)
        indices0.push_back(n0-1);
    if (n1>1)
        indices1.push_back(n1-1);
    if (n2>1)
        indices2.push_back(n2-1);
    for (auto i0:indices0)
        for (auto i1:indices1)
            for (auto i2:indices2)
                positions.push_back((i0*n1+i1)*n2+i2);
    return positions;
}

template<typename T>
inline T sz3_lorenzo_corner_prediction(vector<T> &corners,uint32_t corner_id)
{
    if (corners.size()==4)
    {
        if (corner_id==1 or corner_id==2)
            return corners[0];
        return corners[1]+corners[2]-corners[0];
    }
    if (corner_id==1 or corner_id==2 or corner_id==4)
        return corners[0];
    if (corner_id==3)
        return corners[1]+corners[2]-corners[0];
    if (corner_id==5)
        return corners[1]+corners[4]-corners[0];
    if (corner_id==6)
        return corners[2]+corners[4]-corners[0];
    return corners[3]+corners[5]+corners[6]-corners[1]-corners[2]-corners[4]+corners[0];
}

template<typename T>
void sz3_generate_corner_variant_qbs(vector<T> &data,vector<uint32_t> &positions,uint64_t begin_index,uint64_t end_index,uint32_t element_size,T abs_eb,T neg_2_eb,
    int64_t *independent_previous_qb,int64_t *lorenzo_zero_qb,int64_t *lorenzo_previous_qb,int64_t *mean_zero_qb,int64_t *mean_previous_qb)
{
    uint32_t corner_num=positions.size();
    T two_eb=2*abs_eb;
    vector<T> independent_previous(corner_num,0),lorenzo_zero(corner_num),lorenzo_previous(corner_num);
    T previous_lorenzo_anchor=0,previous_mean=0;
    for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
    {
        T *element=data.data()+element_id*element_size;
        vector<T> truth(corner_num);
        for (uint32_t corner_id=0;corner_id<corner_num;corner_id++)
            truth[corner_id]=element[positions[corner_id]];

        for (uint32_t corner_id=0;corner_id<corner_num;corner_id++)
        {
            T prediction=independent_previous[corner_id];
            int64_t qb=static_cast<int64_t>(round((prediction-truth[corner_id])*neg_2_eb));
            independent_previous_qb[element_id*corner_num+corner_id]=qb;
            independent_previous[corner_id]=prediction-qb*two_eb;

            prediction=corner_id==0?0:sz3_lorenzo_corner_prediction(lorenzo_zero,corner_id);
            qb=static_cast<int64_t>(round((prediction-truth[corner_id])*neg_2_eb));
            lorenzo_zero_qb[element_id*corner_num+corner_id]=qb;
            lorenzo_zero[corner_id]=prediction-qb*two_eb;

            prediction=corner_id==0?previous_lorenzo_anchor:sz3_lorenzo_corner_prediction(lorenzo_previous,corner_id);
            qb=static_cast<int64_t>(round((prediction-truth[corner_id])*neg_2_eb));
            lorenzo_previous_qb[element_id*corner_num+corner_id]=qb;
            lorenzo_previous[corner_id]=prediction-qb*two_eb;
        }
        previous_lorenzo_anchor=lorenzo_previous[0];

        T mean=accumulate(truth.begin(),truth.end(),T(0))/corner_num;
        uint64_t mean_offset=element_id*(corner_num+1);
        int64_t mean_qb=static_cast<int64_t>(round(-mean*neg_2_eb));
        T reconstructed_mean=-mean_qb*two_eb;
        mean_zero_qb[mean_offset]=mean_qb;
        for (uint32_t corner_id=0;corner_id<corner_num;corner_id++)
            mean_zero_qb[mean_offset+1+corner_id]=static_cast<int64_t>(round((reconstructed_mean-truth[corner_id])*neg_2_eb));

        mean_qb=static_cast<int64_t>(round((previous_mean-mean)*neg_2_eb));
        reconstructed_mean=previous_mean-mean_qb*two_eb;
        mean_previous_qb[mean_offset]=mean_qb;
        for (uint32_t corner_id=0;corner_id<corner_num;corner_id++)
            mean_previous_qb[mean_offset+1+corner_id]=static_cast<int64_t>(round((reconstructed_mean-truth[corner_id])*neg_2_eb));
        previous_mean=reconstructed_mean;
    }
}

template<typename Blocks>
pair<uint64_t,uint64_t> sz3_compress_qb_stream(int64_t *qb,uint32_t symbols_per_element,Blocks &blocks,uint32_t element_size,uint32_t max_cpus,bool transpose=false)
{
    vector<uint64_t> huffman_sizes(blocks.size()),zstd_sizes(blocks.size());
    #pragma omp parallel for schedule(static) num_threads(max_cpus)
    for (uint64_t block_id=0;block_id<blocks.size();block_id++)
    {
        auto &block=blocks[block_id];
        uint64_t element_num=block.end_index-block.begin_index;
        uint64_t symbol_num=symbols_per_element*element_num;
        int64_t *block_qb=qb+symbols_per_element*block.begin_index;
        vector<int64_t> transposed_qb;
        if (transpose)
        {
            transposed_qb.resize(symbol_num);
            for (uint32_t symbol_id=0;symbol_id<symbols_per_element;symbol_id++)
                for (uint64_t element_id=0;element_id<element_num;element_id++)
                    transposed_qb[symbol_id*element_num+element_id]=block_qb[element_id*symbols_per_element+symbol_id];
            block_qb=transposed_qb.data();
        }
        vector<unsigned char> huffman_vec(element_size*element_num*sizeof(int64_t));
        unsigned char *huffman=huffman_vec.data();
        mrcz::huffman_encoder_class<int64_t> encoder;
        encoder.preprocess_encode(block_qb,symbol_num);
        encoder.save(huffman);
        encoder.encode(block_qb,symbol_num,huffman);
        huffman_sizes[block_id]=huffman-huffman_vec.data();
        vector<unsigned char> zstd_vec(ZSTD_compressBound(huffman_sizes[block_id]));
        zstd_sizes[block_id]=ZSTD_compress(zstd_vec.data(),zstd_vec.size(),huffman_vec.data(),huffman_sizes[block_id],3);
    }
    return {accumulate(huffman_sizes.begin(),huffman_sizes.end(),uint64_t(0)),accumulate(zstd_sizes.begin(),zstd_sizes.end(),uint64_t(0))};
}

template<typename T>
void sz3_2_compress(mrcz::param_class<T> &params,string field,bool separate_lossless_streams=false)
{
    auto &header=params.header;
    auto &blocks=params.partitions;
    auto &data=params.data;
    auto &x0=params.x[0];
    auto &x1=params.x[1];
    auto &x2=params.x[2];

    uint32_t n0=header.element_shape[0];
    uint32_t n1=header.element_shape[1];
    uint32_t n2=header.element_shape[2];
    uint32_t n12=n1*n2;
    uint32_t n012=n0*n1*n2;
    uint32_t corner_num=1;
    for (uint32_t dim=0;dim<3;dim++)
        if (header.element_shape[dim]>1)
            corner_num*=2;
    T abs_eb=params.abs_eb[field];
    T neg_2_eb=1/(2*abs_eb);
    int64_t *corner_qb=new int64_t[corner_num*params.element_num_to_process];
    int64_t *err_qb=new int64_t[params.element_num_to_process*n012];
    int64_t *no_corner_err_qb=new int64_t[params.element_num_to_process*(n012-corner_num)];
    vector<uint32_t> corner_positions=sz3_corner_positions(n0,n1,n2);
    vector<int64_t> independent_previous_qb,lorenzo_zero_qb,lorenzo_previous_qb,mean_zero_qb,mean_previous_qb;
    if (separate_lossless_streams)
    {
        uint64_t corner_count=corner_num*params.element_num_to_process;
        uint64_t mean_count=(corner_num+1)*params.element_num_to_process;
        independent_previous_qb.resize(corner_count);
        lorenzo_zero_qb.resize(corner_count);
        lorenzo_previous_qb.resize(corner_count);
        mean_zero_qb.resize(mean_count);
        mean_previous_qb.resize(mean_count);
    }
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
        sz3_generate_weight(left_indices0,right_indices0,mid_indices0,compressed_indices0,left_weight0,right_weight0,x0);
    }
    if (n1==1)
        compressed_indices1.push_back({0});
    else
    {
        compressed_indices1.push_back({0,n1-1});
        sz3_generate_weight(left_indices1,right_indices1,mid_indices1,compressed_indices1,left_weight1,right_weight1,x1);
    }
    if (n2==1)
        compressed_indices2.push_back({0});
    else
    {
        compressed_indices2.push_back({0,n2-1});
        sz3_generate_weight(left_indices2,right_indices2,mid_indices2,compressed_indices2,left_weight2,right_weight2,x2);
    }
    uint32_t max_level=max(max(left_indices0.size(),left_indices1.size()),left_indices2.size());
    vector<uint64_t> corner_huffman_sizes(blocks.size()),corner_zstd_sizes(blocks.size());
    vector<uint64_t> residual_huffman_sizes(blocks.size()),residual_zstd_sizes(blocks.size());

    double start_time=omp_get_wtime();
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
    for (uint64_t block_id=0;block_id<blocks.size();block_id++) //这个block_id是假的，指的是在需要处理的blocks中的编号，而不是header文件中的编号
    {
        auto &block=blocks[block_id];
        uint64_t begin_index=block.begin_index;
        uint64_t end_index=block.end_index;
        if (separate_lossless_streams)
            sz3_generate_corner_variant_qbs(data,corner_positions,begin_index,end_index,n012,abs_eb,neg_2_eb,independent_previous_qb.data(),lorenzo_zero_qb.data(),lorenzo_previous_qb.data(),mean_zero_qb.data(),mean_previous_qb.data());
        // 先量化位于每个element位于角上的点，覆盖原数据
        // 之后可能要改成压缩率更高的方法
        params.quantize_and_overwrite_corner_numbers(corner_qb,begin_index,end_index,n0,n1,n2,n12,n012,abs_eb,neg_2_eb);
        // 然后量化剩余的数据
        for (uint64_t i=begin_index;i<end_index;i++)
        {
            T *element=data.data()+i*n012;
            int64_t *element_err_qb=err_qb+i*n012;
            int64_t *element_no_corner_err_qb=no_corner_err_qb+i*(n012-corner_num);
            uint32_t level0=0,level1=0,level2=0;
            for (uint32_t level=0;level<max_level;level++)
            {
                if (n0>1 and level0<left_indices0.size())
                {
                    for (uint32_t i0=0;i0<left_indices0[level0].size();i0++)
                        for (auto i1:compressed_indices1[level1])
                            for (auto i2:compressed_indices2[level2])
                            {
                                uint32_t pos_base=i1*n2+i2;
                                uint32_t mid_pos=pos_base+mid_indices0[level0][i0]*n12;
                                T left_value=element[pos_base+left_indices0[level0][i0]*n12];
                                T right_value=element[pos_base+right_indices0[level0][i0]*n12];
                                T pred_value=left_weight0[level0][i0]*left_value+right_weight0[level0][i0]*right_value;
                                params.quantize_and_overwrite_err(pred_value,element[mid_pos],element_err_qb[mid_pos],abs_eb,neg_2_eb);
                            }
                    level0++;
                }
                if (n1>1 and level1<left_indices1.size())
                {
                    for (auto i0:compressed_indices0[level0])
                        for (uint32_t i1=0;i1<left_indices1[level1].size();i1++)
                            for (auto i2:compressed_indices2[level2])
                            {
                                uint32_t pos_base=i0*n12+i2;
                                uint32_t mid_pos=pos_base+mid_indices1[level1][i1]*n2;
                                T left_value=element[pos_base+left_indices1[level1][i1]*n2];
                                T right_value=element[pos_base+right_indices1[level1][i1]*n2];
                                T pred_value=left_weight1[level1][i1]*left_value+right_weight1[level1][i1]*right_value;
                                params.quantize_and_overwrite_err(pred_value,element[mid_pos],element_err_qb[mid_pos],abs_eb,neg_2_eb);
                            }
                    level1++;
                }
                if (n2>1 and level2<left_indices2.size())
                {
                    for (auto i0:compressed_indices0[level0])
                        for (auto i1:compressed_indices1[level1])
                            for (uint32_t i2=0;i2<left_indices2[level2].size();i2++)
                            {
                                uint32_t pos_base=i0*n12+i1*n2;
                                uint32_t mid_pos=pos_base+mid_indices2[level2][i2];
                                T left_value=element[pos_base+left_indices2[level2][i2]];
                                T right_value=element[pos_base+right_indices2[level2][i2]];
                                T pred_value=left_weight2[level2][i2]*left_value+right_weight2[level2][i2]*right_value;
                                params.quantize_and_overwrite_err(pred_value,element[mid_pos],element_err_qb[mid_pos],abs_eb,neg_2_eb);
                            }
                    level2++;
                }
            }
            params.convert_err_eb_to_no_corner_err_eb(element_no_corner_err_qb,element_err_qb,n0,n1,n2,n12,n012);
        }
        if (separate_lossless_streams)
        {
            uint64_t element_num=end_index-begin_index;
            vector<unsigned char> corner_huffman_vec(n012*element_num*sizeof(int64_t));
            unsigned char *corner_huffman=corner_huffman_vec.data();
            mrcz::huffman_encoder_class<int64_t> corner_encoder;
            corner_encoder.preprocess_encode(corner_qb+corner_num*begin_index,corner_num*element_num);
            corner_encoder.save(corner_huffman);
            corner_encoder.encode(corner_qb+corner_num*begin_index,corner_num*element_num,corner_huffman);
            corner_huffman_sizes[block_id]=corner_huffman-corner_huffman_vec.data();
            vector<unsigned char> corner_zstd_vec(ZSTD_compressBound(corner_huffman_sizes[block_id]));
            corner_zstd_sizes[block_id]=ZSTD_compress(corner_zstd_vec.data(),corner_zstd_vec.size(),corner_huffman_vec.data(),corner_huffman_sizes[block_id],3);

            vector<unsigned char> residual_huffman_vec(n012*element_num*sizeof(int64_t));
            unsigned char *residual_huffman=residual_huffman_vec.data();
            mrcz::huffman_encoder_class<int64_t> residual_encoder;
            residual_encoder.preprocess_encode(no_corner_err_qb+begin_index*(n012-corner_num),(n012-corner_num)*element_num);
            residual_encoder.save(residual_huffman);
            residual_encoder.encode(no_corner_err_qb+begin_index*(n012-corner_num),(n012-corner_num)*element_num,residual_huffman);
            residual_huffman_sizes[block_id]=residual_huffman-residual_huffman_vec.data();
            vector<unsigned char> residual_zstd_vec(ZSTD_compressBound(residual_huffman_sizes[block_id]));
            residual_zstd_sizes[block_id]=ZSTD_compress(residual_zstd_vec.data(),residual_zstd_vec.size(),residual_huffman_vec.data(),residual_huffman_sizes[block_id],3);
            continue;
        }
        // 把两种qb无损压缩到一起
        vector<unsigned char> huffman_bs_vec(n012*(end_index-begin_index)*sizeof(int64_t));
        unsigned char *huffman_bs=huffman_bs_vec.data();
        mrcz::huffman_encoder_class<int64_t> corner_encoder;
        corner_encoder.preprocess_encode(corner_qb+corner_num*begin_index,corner_num*(end_index-begin_index));
        corner_encoder.save(huffman_bs);
        corner_encoder.encode(corner_qb+corner_num*begin_index,corner_num*(end_index-begin_index),huffman_bs);
        mrcz::huffman_encoder_class<int64_t> err_encoder;
        err_encoder.preprocess_encode(no_corner_err_qb+begin_index*(n012-corner_num),(n012-corner_num)*(end_index-begin_index));
        err_encoder.save(huffman_bs);
        err_encoder.encode(no_corner_err_qb+begin_index*(n012-corner_num),(n012-corner_num)*(end_index-begin_index),huffman_bs);
        uint64_t huffman_bs_len=huffman_bs-huffman_bs_vec.data();
        uint64_t max_zstd_bs_len=ZSTD_compressBound(huffman_bs_len);
        block.zstd_bs.resize(max_zstd_bs_len);
        uint64_t zstd_bs_len=ZSTD_compress(block.zstd_bs.data(),max_zstd_bs_len,huffman_bs_vec.data(),huffman_bs_len,3);
        block.zstd_bs.resize(zstd_bs_len);
    }
    double end_time=omp_get_wtime();
    printf("Compression time= %f seconds\n",end_time-start_time);
    printf("\033[91mField= %s, compression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),params.element_num_to_process*n012*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
    if (separate_lossless_streams)
    {
        uint64_t corner_huffman_size=accumulate(corner_huffman_sizes.begin(),corner_huffman_sizes.end(),uint64_t(0));
        uint64_t corner_zstd_size=accumulate(corner_zstd_sizes.begin(),corner_zstd_sizes.end(),uint64_t(0));
        uint64_t residual_huffman_size=accumulate(residual_huffman_sizes.begin(),residual_huffman_sizes.end(),uint64_t(0));
        uint64_t residual_zstd_size=accumulate(residual_zstd_sizes.begin(),residual_zstd_sizes.end(),uint64_t(0));
        double value_num=params.element_num_to_process*n012;
        uint64_t total_zstd_size=corner_zstd_size+residual_zstd_size;
        uint64_t corner_count=corner_num*params.element_num_to_process;
        uint64_t residual_count=(n012-corner_num)*params.element_num_to_process;
        uint64_t mean_count=(corner_num+1)*params.element_num_to_process;
        auto residual_statistics=sz3_quantization_statistics(no_corner_err_qb,residual_count);
        auto independent_previous_sizes=sz3_compress_qb_stream(independent_previous_qb.data(),corner_num,blocks,n012,params.max_cpus);
        auto lorenzo_zero_sizes=sz3_compress_qb_stream(lorenzo_zero_qb.data(),corner_num,blocks,n012,params.max_cpus);
        auto lorenzo_previous_sizes=sz3_compress_qb_stream(lorenzo_previous_qb.data(),corner_num,blocks,n012,params.max_cpus);
        auto mean_zero_sizes=sz3_compress_qb_stream(mean_zero_qb.data(),corner_num+1,blocks,n012,params.max_cpus);
        auto mean_previous_sizes=sz3_compress_qb_stream(mean_previous_qb.data(),corner_num+1,blocks,n012,params.max_cpus);
        auto independent_zero_transposed_sizes=sz3_compress_qb_stream(corner_qb,corner_num,blocks,n012,params.max_cpus,true);
        auto independent_previous_transposed_sizes=sz3_compress_qb_stream(independent_previous_qb.data(),corner_num,blocks,n012,params.max_cpus,true);
        auto lorenzo_zero_transposed_sizes=sz3_compress_qb_stream(lorenzo_zero_qb.data(),corner_num,blocks,n012,params.max_cpus,true);
        auto lorenzo_previous_transposed_sizes=sz3_compress_qb_stream(lorenzo_previous_qb.data(),corner_num,blocks,n012,params.max_cpus,true);
        auto mean_zero_transposed_sizes=sz3_compress_qb_stream(mean_zero_qb.data(),corner_num+1,blocks,n012,params.max_cpus,true);
        auto mean_previous_transposed_sizes=sz3_compress_qb_stream(mean_previous_qb.data(),corner_num+1,blocks,n012,params.max_cpus,true);
        auto print_corner_variant=[&](string name,int64_t *qb,uint64_t symbol_count,uint64_t huffman_size,uint64_t zstd_size)
        {
            auto statistics=sz3_quantization_statistics(qb,symbol_count);
            printf("Corner variant %-22s Huffman=%lu bytes, Zstd=%lu bytes, bits/original-corner=%f, bits/stored-symbol=%f, zero=%f%%, entropy=%f bits/symbol\n",name.c_str(),huffman_size,zstd_size,zstd_size*8.0/corner_count,zstd_size*8.0/symbol_count,statistics.first*100,statistics.second);
        };
        printf("\033[91m");
        print_corner_variant("independent+zero",corner_qb,corner_count,corner_huffman_size,corner_zstd_size);
        print_corner_variant("independent+previous",independent_previous_qb.data(),corner_count,independent_previous_sizes.first,independent_previous_sizes.second);
        print_corner_variant("Lorenzo+zero",lorenzo_zero_qb.data(),corner_count,lorenzo_zero_sizes.first,lorenzo_zero_sizes.second);
        print_corner_variant("Lorenzo+previous",lorenzo_previous_qb.data(),corner_count,lorenzo_previous_sizes.first,lorenzo_previous_sizes.second);
        print_corner_variant("mean+zero",mean_zero_qb.data(),mean_count,mean_zero_sizes.first,mean_zero_sizes.second);
        print_corner_variant("mean+previous",mean_previous_qb.data(),mean_count,mean_previous_sizes.first,mean_previous_sizes.second);
        print_corner_variant("independent+zero+T",corner_qb,corner_count,independent_zero_transposed_sizes.first,independent_zero_transposed_sizes.second);
        print_corner_variant("independent+previous+T",independent_previous_qb.data(),corner_count,independent_previous_transposed_sizes.first,independent_previous_transposed_sizes.second);
        print_corner_variant("Lorenzo+zero+T",lorenzo_zero_qb.data(),corner_count,lorenzo_zero_transposed_sizes.first,lorenzo_zero_transposed_sizes.second);
        print_corner_variant("Lorenzo+previous+T",lorenzo_previous_qb.data(),corner_count,lorenzo_previous_transposed_sizes.first,lorenzo_previous_transposed_sizes.second);
        print_corner_variant("mean+zero+T",mean_zero_qb.data(),mean_count,mean_zero_transposed_sizes.first,mean_zero_transposed_sizes.second);
        print_corner_variant("mean+previous+T",mean_previous_qb.data(),mean_count,mean_previous_transposed_sizes.first,mean_previous_transposed_sizes.second);
        printf("Current corner share=%f%%\n",corner_zstd_size*100.0/total_zstd_size);
        printf("Residual stream: Huffman=%lu bytes, Zstd=%lu bytes, BR=%f bits/value, bits/non-corner=%f, share=%f%%, zero=%f%%, entropy=%f bits/symbol\n",residual_huffman_size,residual_zstd_size,residual_zstd_size*8/value_num,residual_zstd_size*8.0/residual_count,residual_zstd_size*100.0/total_zstd_size,residual_statistics.first*100,residual_statistics.second);
        printf("Separate streams total: Zstd=%lu bytes, BR=%f bits/value\n\033[0m",total_zstd_size,total_zstd_size*8/value_num);
    }
    delete[] corner_qb;
    delete[] err_qb;
    delete[] no_corner_err_qb;
}

template<typename T>
void sz3_3d_gll_experimental_compress(mrcz::param_class<T> &params,string field)
{
    sz3_2_compress(params,field,true);
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

