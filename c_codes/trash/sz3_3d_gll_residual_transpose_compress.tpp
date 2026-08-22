// XXX 这是全部非角点prediction residual跨element转置实验的压缩端快照。
// XXX 完整测试显示主要数据的BR和速度均变差，因此不参与构建。

#pragma once


#include <set>
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

template<typename T>
void sz3_2_compress(mrcz::param_class<T> &params,string field,bool transpose_residual=false)
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
    uint32_t residual_num=n012-corner_num;
    T abs_eb=params.abs_eb[field];
    T neg_2_eb=1/(2*abs_eb);
    int64_t *corner_qb=new int64_t[corner_num*params.element_num_to_process];
    int64_t *err_qb=new int64_t[params.element_num_to_process*n012];
    int64_t *no_corner_err_qb=new int64_t[params.element_num_to_process*residual_num];
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
    double start_time=omp_get_wtime();
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
    for (uint64_t block_id=0;block_id<blocks.size();block_id++) //这个block_id是假的，指的是在需要处理的blocks中的编号，而不是header文件中的编号
    {
        auto &block=blocks[block_id];
        uint64_t begin_index=block.begin_index;
        uint64_t end_index=block.end_index;
        uint64_t element_num=end_index-begin_index;
        int64_t *block_no_corner_err_qb=no_corner_err_qb+begin_index*residual_num;
        vector<int64_t> element_no_corner_err_qb_buffer;
        if (transpose_residual)
            element_no_corner_err_qb_buffer.resize(residual_num);
        // 先量化位于每个element位于角上的点，覆盖原数据
        // 之后可能要改成压缩率更高的方法
        params.quantize_and_overwrite_corner_numbers(corner_qb,begin_index,end_index,n0,n1,n2,n12,n012,abs_eb,neg_2_eb);
        // 然后量化剩余的数据
        for (uint64_t i=begin_index;i<end_index;i++)
        {
            T *element=data.data()+i*n012;
            int64_t *element_err_qb=err_qb+i*n012;
            int64_t *element_no_corner_err_qb=transpose_residual?element_no_corner_err_qb_buffer.data():no_corner_err_qb+i*residual_num;
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
            if (transpose_residual)
            {
                uint64_t local_element_id=i-begin_index;
                for (uint32_t residual_id=0;residual_id<residual_num;residual_id++)
                    block_no_corner_err_qb[residual_id*element_num+local_element_id]=element_no_corner_err_qb[residual_id];
            }
        }
        int64_t *block_corner_qb=corner_qb+corner_num*begin_index;
        // 把两种qb无损压缩到一起
        vector<unsigned char> huffman_bs_vec(n012*element_num*sizeof(int64_t));
        unsigned char *huffman_bs=huffman_bs_vec.data();
        mrcz::huffman_encoder_class<int64_t> corner_encoder;
        corner_encoder.preprocess_encode(block_corner_qb,corner_num*element_num);
        corner_encoder.save(huffman_bs);
        corner_encoder.encode(block_corner_qb,corner_num*element_num,huffman_bs);
        mrcz::huffman_encoder_class<int64_t> err_encoder;
        err_encoder.preprocess_encode(block_no_corner_err_qb,residual_num*element_num);
        err_encoder.save(huffman_bs);
        err_encoder.encode(block_no_corner_err_qb,residual_num*element_num,huffman_bs);
        uint64_t huffman_bs_len=huffman_bs-huffman_bs_vec.data();
        uint64_t max_zstd_bs_len=ZSTD_compressBound(huffman_bs_len);
        block.zstd_bs.resize(max_zstd_bs_len);
        uint64_t zstd_bs_len=ZSTD_compress(block.zstd_bs.data(),max_zstd_bs_len,huffman_bs_vec.data(),huffman_bs_len,3);
        block.zstd_bs.resize(zstd_bs_len);
    }
    double end_time=omp_get_wtime();
    printf("Compression time= %f seconds\n",end_time-start_time);
    printf("\033[91mField= %s, compression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),params.element_num_to_process*n012*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
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

