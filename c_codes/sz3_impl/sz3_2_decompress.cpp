#pragma once

#include "../mrcz_core/mrcz_param.hpp"
#include "sz3_2_compress.cpp"

template<typename T>
void sz3_2_decompress(mrcz::param_class<T> &params,string field)
{
	auto &header=params.header;
    auto &blocks=params.blocks;
    auto &decompressed_data=params.decompressed_data;
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
    int64_t *corner_qb=new int64_t[8*params.element_num_to_process];
    int64_t *err_qb=new int64_t[params.element_num_to_process*n012];
    int64_t *no_corner_err_qb=new int64_t[params.element_num_to_process*(n012-corner_num)];
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
    for (uint64_t block_id=0;block_id<blocks.size();block_id++)
    {
        auto &block=blocks[block_id];
        uint64_t begin_index=block.begin_index;
        uint64_t end_index=block.end_index;
        // 先无损解压
        uint64_t huffman_bs_len=ZSTD_getFrameContentSize(block.zstd_bs.data(),block.zstd_bs.size());
        vector<unsigned char> huffman_bs_vec(huffman_bs_len);
        if (ZSTD_decompress(huffman_bs_vec.data(),huffman_bs_len,block.zstd_bs.data(),block.zstd_bs.size())!=huffman_bs_len)
            throw runtime_error("ZSTD_decompress error in decompress()");
        const unsigned char *huffman_bs=huffman_bs_vec.data();
        mrcz::huffman_encoder_class<int64_t> corner_encoder;
        corner_encoder.load(huffman_bs);
        corner_encoder.decode(huffman_bs,corner_qb+8*begin_index,8*(end_index-begin_index));
        mrcz::huffman_encoder_class<int64_t> err_encoder;
        err_encoder.load(huffman_bs);
        err_encoder.decode(huffman_bs,no_corner_err_qb+(n012-corner_num)*begin_index,(n012-corner_num)*(end_index-begin_index));
        // 重构每个element角上的数据
        params.reconstruct_corner_number(corner_qb,begin_index,end_index,n2,n12,n012,abs_eb);
        //把no_corner_err_qb转换成err_qb
        for (uint64_t i=begin_index;i<end_index;i++)
        {
            int64_t *element_err_qb=err_qb+i*n012;
            int64_t *element_no_corner_err_qb=no_corner_err_qb+i*(n012-corner_num);
            if (n2>1)
            {
                memcpy(element_err_qb+1,element_no_corner_err_qb,sizeof(int64_t)*(n2-2));
                memcpy(element_err_qb+n2,element_no_corner_err_qb+n2-2,sizeof(int64_t)*(n12-2*n2));
                memcpy(element_err_qb+n12-n2+1,element_no_corner_err_qb+n12-n2-2,sizeof(int64_t)*(n2-2));
                memcpy(element_err_qb+n12,element_no_corner_err_qb+n12-4,sizeof(int64_t)*(n012-2*n12));
                memcpy(element_err_qb+n012-n12+1,element_no_corner_err_qb+n012-n12-4,sizeof(int64_t)*(n2-2));
                memcpy(element_err_qb+n012-n12+n2,element_no_corner_err_qb+n012-n12+n2-6,sizeof(int64_t)*(n12-2*n2));
                memcpy(element_err_qb+n012-n2+1,element_no_corner_err_qb+n012-n2-6,sizeof(int64_t)*(n2-2));
            }
            else
            {
                memcpy(element_err_qb+n2,element_no_corner_err_qb+n2-1,sizeof(int64_t)*(n12-2*n2));
                memcpy(element_err_qb+n12,element_no_corner_err_qb+n12-2,sizeof(int64_t)*(n012-2*n12));
                memcpy(element_err_qb+n012-n12+n2,element_no_corner_err_qb+n012-n12+n2-3,sizeof(int64_t)*(n12-2*n2));
            }
            // for (uint32_t i=0;i<n012;i++)
            //     printf("%ld ",element_err_qb[i]);
            // printf("\n");

            T *element=decompressed_data.data()+i*n012;
            uint32_t level0=0,level1=0,level2=0;
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
                                element[mid_pos]=pred_value-2*abs_eb*element_err_qb[mid_pos];
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
                                element[mid_pos]=pred_value-2*abs_eb*element_err_qb[mid_pos];
                                // printf("mid_pos=%u, pred_value=%f, err_qb=%ld, decompressed_value=%f\n",mid_pos,pred_value,element_err_qb[mid_pos],element[mid_pos]);
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
                                element[mid_pos]=pred_value-2*abs_eb*element_err_qb[mid_pos];
                                // printf("mid_pos=%u, pred_value=%f, err_qb=%ld, decompressed_value=%f\n",mid_pos,pred_value,element_err_qb[mid_pos],element[mid_pos]);
                            }
                    level0++;
                }
                // printf("ccc\n");
            }
        }
    }
    double end_time=omp_get_wtime();
    delete[] corner_qb;
    delete[] err_qb;
    delete[] no_corner_err_qb;
    printf("Decompression time= %f seconds\n",end_time-start_time);
    printf("\033[91mBandwidth per cpu= %f MB/s\n\033[0m",params.element_num_to_process*n012*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
}