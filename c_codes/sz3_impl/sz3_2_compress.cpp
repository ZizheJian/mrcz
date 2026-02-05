#pragma once

#include <set>
#include "../mrcz_core/mrcz_param.hpp"
#include "../mrcz_core/huffman.hpp"
#include "../mrcz_core/quantize.tpp"
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
        for (uint32_t i=0;i<left_indices.back().size();i++)
        {
            uint32_t left_index=left_indices.back()[i],right_index=right_indices.back()[i],mid_index=left_index;
            T left_pos=x[left_index],right_pos=x[right_index],ideal_mid_pos=(left_pos+right_pos)/2;
            for (uint32_t j=left_index+1;j<right_index;j++)
                if (fabs(x[j]-ideal_mid_pos)<fabs(x[mid_index]-ideal_mid_pos))
                    mid_index=j;
            mid_indices.back().push_back(mid_index);
            T mid_pos=x[mid_index];
            T left_dist=mid_pos-left_pos;
            T right_dist=right_pos-mid_pos;
            left_weight.back().push_back(right_dist/(left_dist+right_dist));
            right_weight.back().push_back(left_dist/(left_dist+right_dist));
        }
        vector<uint32_t> new_left_indices,new_right_indices;
        for (uint32_t i=0;i<left_indices.back().size();i++)
        {
            uint32_t left_index=left_indices.back()[i],right_index=right_indices.back()[i],mid_index=mid_indices.back()[i];
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
void sz3_2_compress(mrcz::param_class<T> &params,string field)
{
    auto &header=params.header;
    auto &blocks=params.blocks;
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
    int64_t *corner_qb=new int64_t[8*params.element_num_to_process];
    int64_t *err_qb=new int64_t[params.element_num_to_process*n012];
    int64_t *no_corner_err_qb=new int64_t[params.element_num_to_process*(n012-corner_num)];
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
    // for (int level=0;level<left_indices0.size();level++)
    //     for (int i=0;i<left_indices0[level].size();i++)
    //         printf("level0=%d, i=%d, left_index=%u, right_index=%u, mid_index=%u, left_weight=%f, right_weight=%f\n",level,i,
    //             left_indices0[level][i],right_indices0[level][i],mid_indices0[level][i],
    //             left_weight0[level][i],right_weight0[level][i]);
    // for (int level=0;level<left_indices1.size();level++)
    //     for (int i=0;i<left_indices1[level].size();i++)
    //         printf("level1=%d, i=%d, left_index=%u, right_index=%u, mid_index=%u, left_weight=%f, right_weight=%f\n",level,i,
    //             left_indices1[level][i],right_indices1[level][i],mid_indices1[level][i],
    //             left_weight1[level][i],right_weight1[level][i]);
    uint32_t max_level=max(max(left_indices0.size(),left_indices1.size()),left_indices2.size());

    double start_time=omp_get_wtime();
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
    for (uint64_t block_id=0;block_id<blocks.size();block_id++) //这个block_id是假的，指的是在需要处理的blocks中的编号，而不是header文件中的编号
    {
        auto &block=blocks[block_id];
        uint64_t begin_index=block.begin_index;
        uint64_t end_index=block.end_index;
        // 先量化位于每个element位于角上的点，覆盖原数据
        // 之后可能要改成压缩率更高的方法
        params.quantize_and_overwrite_corner_numbers(corner_qb,begin_index,end_index,n2,n12,n012,abs_eb,neg_2_eb);
        // 然后量化剩余的数据
        for (uint64_t i=begin_index;i<end_index;i++)
        {
            T *element=data.data()+i*n012;
            int64_t *element_err_qb=err_qb+i*n012;
            int64_t *element_no_corner_err_qb=no_corner_err_qb+i*(n012-corner_num);
            uint32_t level0=0,level1=0,level2=0;
            // printf("max_level=%u\n",max_level);
            for (uint32_t level=0;level<max_level;level++)
            {
                // printf("level=%u, level0=%u, level1=%u, level2=%u\n",level,level0,level1,level2);
                if (n2>1 && level2<left_indices2.size())
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
                // printf("aaa\n");
                if (n1>1 && level1<left_indices1.size())
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
                                // printf("mid_pos=%u, pred_value=%f, original_value=%f\n",mid_pos,pred_value,element[mid_pos]);
                                params.quantize_and_overwrite_err(pred_value,element[mid_pos],element_err_qb[mid_pos],abs_eb,neg_2_eb);
                                // printf("element_err_qb=%ld\n",element_err_qb[mid_pos]);
                            }
                    level1++;
                }
                // printf("bbb\n");
                if (n0>1 && level0<left_indices0.size())
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
                                // printf("mid_pos=%u, pred_value=%f, original_value=%f\n",mid_pos,pred_value,element[mid_pos]);
                                params.quantize_and_overwrite_err(pred_value,element[mid_pos],element_err_qb[mid_pos],abs_eb,neg_2_eb);
                                // printf("element_err_qb=%ld\n",element_err_qb[mid_pos]);
                            }
                    level0++;
                }
                // printf("ccc\n");
            }
            // for (uint32_t i=0;i<n012;i++)
            //     printf("%ld ",element_err_qb[i]);
            // printf("\n");
            params.convert_err_eb_to_no_corner_err_eb(element_no_corner_err_qb,element_err_qb,n2,n12,n012);
        }
        // 把两种qb无损压缩到一起
        vector<unsigned char> huffman_bs_vec((2*n012+8)*(end_index-begin_index)*sizeof(int64_t));
        unsigned char *huffman_bs=huffman_bs_vec.data();
        mrcz::huffman_encoder_class<int64_t> corner_encoder;
        corner_encoder.preprocess_encode(corner_qb+8*begin_index,8*(end_index-begin_index));
        corner_encoder.save(huffman_bs);
        corner_encoder.encode(corner_qb+8*begin_index,8*(end_index-begin_index),huffman_bs);
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
    delete[] corner_qb;
    delete[] err_qb;
    delete[] no_corner_err_qb;
    printf("Compression time= %f seconds\n",end_time-start_time);
    printf("\033[91mBandwidth per cpu= %f MB/s\n\033[0m",params.element_num_to_process*n012*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
}