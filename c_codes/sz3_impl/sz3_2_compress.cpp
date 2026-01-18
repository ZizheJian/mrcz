#pragma once

#include <set>
#include "../jzzz_core/jzzz_param.hpp"
#include "../jzzz_core/huffman.hpp"
#include "../jzzz_core/quantize.tpp"
#include "../third_party/zstd/lib/zstd.h"


template<typename T>
void sz3_2_compress(jzzz::param_class<T> &params)
{
    auto &header=params.header;
    auto &blocks=params.blocks;
    auto &data=params.data;
    auto &zstd_bs=params.zstd_bs;
    auto &abs_eb=params.abs_eb;
    auto &x=params.x;

	const char *slurm_cpus=getenv("SLURM_CPUS_PER_TASK");
    int available_cpus=0;
    if (slurm_cpus!=nullptr)
        available_cpus=atoi(slurm_cpus);
    else
        throw runtime_error("没有获取到每个task可用的cpu数量");
    printf("Available cpus per task= %d\n",available_cpus);
    uint32_t n=header.element_shape[0];
    uint32_t n2=n*n;
    uint32_t n3=n2*n;
    T neg_2_eb=1/(2*abs_eb);
    int64_t *corner_qb=new int64_t[8*params.element_num_to_process];
    int64_t *err_qb=new int64_t[params.element_num_to_process*header.element_size];
    int64_t *no_corner_err_qb=new int64_t[params.element_num_to_process*(header.element_size-8)];
    zstd_bs.resize(blocks.size());
    vector<vector<uint32_t>> left_indices,right_indices,mid_indices,compressed_indices;
    vector<vector<T>> left_weight,right_weight;
    left_indices.push_back({0});
    right_indices.push_back({n-1});
    compressed_indices.push_back({0,n-1});
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
    double start_time=omp_get_wtime();
    #pragma omp parallel for schedule(static) num_threads(available_cpus)
    for (uint64_t fake_block_id=0;fake_block_id<params.blocks_to_process.size();fake_block_id++) //这个block_id是假的，指的是在需要处理的blocks中的编号，而不是header文件中的编号
    {
        uint64_t block_id=params.blocks_to_process[fake_block_id];
        // printf("CPU %d compressing block %ld\n",omp_get_thread_num(),block_id);
        uint64_t begin_index=blocks[block_id].begin_index;
        uint64_t end_index=blocks[block_id].end_index;
        // 先量化位于每个element位于角上的点，覆盖原数据
        // 之后可能要改成压缩率更高的方法
        vector<T> corner_numbers(8*(end_index-begin_index));
        for (uint64_t i=begin_index;i<end_index;i++)
        {
            T *element=data.data()+i*header.element_size;
            int64_t *element_corner_qb=corner_qb+8*i;
            params.quantize_and_overwrite(element[0],element_corner_qb[0],abs_eb,neg_2_eb);
            params.quantize_and_overwrite(element[n-1],element_corner_qb[1],abs_eb,neg_2_eb);
            params.quantize_and_overwrite(element[n2-n],element_corner_qb[2],abs_eb,neg_2_eb);
            params.quantize_and_overwrite(element[n2-1],element_corner_qb[3],abs_eb,neg_2_eb);
            params.quantize_and_overwrite(element[n3-n2],element_corner_qb[4],abs_eb,neg_2_eb);
            params.quantize_and_overwrite(element[n3-n2+n-1],element_corner_qb[5],abs_eb,neg_2_eb);
            params.quantize_and_overwrite(element[n3-n],element_corner_qb[6],abs_eb,neg_2_eb);
            params.quantize_and_overwrite(element[n3-1],element_corner_qb[7],abs_eb,neg_2_eb);
        }
        // 然后量化剩余的数据
        for (uint64_t i=begin_index;i<end_index;i++)
        {
            T *element=data.data()+i*header.element_size;
            int64_t *element_err_qb=err_qb+i*header.element_size;
            for (uint32_t level=0;level<left_indices.size();level++)
            {
                for (auto i0:compressed_indices[level])
                    for (auto i1:compressed_indices[level])
                        for (uint32_t i2=0;i2<left_indices[level].size();i2++)
                        {
                            uint32_t pos_base=i0*n2+i1*n;
                            uint32_t mid_pos=pos_base+mid_indices[level][i2];
                            T left_value=element[pos_base+left_indices[level][i2]];
                            T right_value=element[pos_base+right_indices[level][i2]];
                            T pred_value=left_weight[level][i2]*left_value+right_weight[level][i2]*right_value;
                            params.quantize_and_overwrite_err(pred_value,element[mid_pos],element_err_qb[mid_pos],abs_eb,neg_2_eb);
                        }
                for (auto i0:compressed_indices[level])
                    for (uint32_t i1=0;i1<left_indices[level].size();i1++)
                        for (auto i2:compressed_indices[level+1])
                        {
                            uint32_t pos_base=i0*n2+i2;
                            uint32_t mid_pos=pos_base+mid_indices[level][i1]*n;
                            T left_value=element[pos_base+left_indices[level][i1]*n];
                            T right_value=element[pos_base+right_indices[level][i1]*n];
                            T pred_value=left_weight[level][i1]*left_value+right_weight[level][i1]*right_value;
                            params.quantize_and_overwrite_err(pred_value,element[mid_pos],element_err_qb[mid_pos],abs_eb,neg_2_eb);
                        }
                for (uint32_t i0=0;i0<left_indices[level].size();i0++)
                    for (auto i1:compressed_indices[level+1])
                        for (auto i2:compressed_indices[level+1])
                        {
                            uint32_t pos_base=i1*n+i2;
                            uint32_t mid_pos=pos_base+mid_indices[level][i0]*n2;
                            T left_value=element[pos_base+left_indices[level][i0]*n2];
                            T right_value=element[pos_base+right_indices[level][i0]*n2];
                            T pred_value=left_weight[level][i0]*left_value+right_weight[level][i0]*right_value;
                            params.quantize_and_overwrite_err(pred_value,element[mid_pos],element_err_qb[mid_pos],abs_eb,neg_2_eb);
                        }
            }
        }
        //去掉err_qb中位于角上的qb
        for (uint64_t i=begin_index;i<end_index;i++)
        {
            int64_t *element_err_qb=err_qb+i*header.element_size;
            int64_t *element_no_corner_err_qb=no_corner_err_qb+i*(header.element_size-8);
            memcpy(element_no_corner_err_qb,element_err_qb+1,sizeof(int64_t)*(n-2));
            memcpy(element_no_corner_err_qb+n-2,element_err_qb+n,sizeof(int64_t)*(n2-2*n));
            memcpy(element_no_corner_err_qb+n2-n-2,element_err_qb+n2-n+1,sizeof(int64_t)*(n-2));
            memcpy(element_no_corner_err_qb+n2-4,element_err_qb+n2,sizeof(int64_t)*(n3-2*n2));
            memcpy(element_no_corner_err_qb+n3-n2-4,element_err_qb+n3-n2+1,sizeof(int64_t)*(n-2));
            memcpy(element_no_corner_err_qb+n3-n2+n-6,element_err_qb+n3-n2+n,sizeof(int64_t)*(n2-2*n));
            memcpy(element_no_corner_err_qb+n3-n-6,element_err_qb+n3-n+1,sizeof(int64_t)*(n-2));
        }
        // 把两种qb无损压缩到一起
        vector<unsigned char> huffman_bs_vec((2*header.element_size+8)*(end_index-begin_index)*sizeof(int64_t));
        unsigned char *huffman_bs=huffman_bs_vec.data();
        jzzz::huffman_encoder_class<int64_t> corner_encoder;
        corner_encoder.preprocess_encode(corner_qb+8*begin_index,8*(end_index-begin_index));
        corner_encoder.save(huffman_bs);
        corner_encoder.encode(corner_qb+8*begin_index,8*(end_index-begin_index),huffman_bs);
        jzzz::huffman_encoder_class<int64_t> err_encoder;
        err_encoder.preprocess_encode(no_corner_err_qb+begin_index*(header.element_size-8),(header.element_size-8)*(end_index-begin_index));
        err_encoder.save(huffman_bs);
        err_encoder.encode(no_corner_err_qb+begin_index*(header.element_size-8),(header.element_size-8)*(end_index-begin_index),huffman_bs);
        uint64_t huffman_bs_len=huffman_bs-huffman_bs_vec.data();
        uint64_t max_zstd_bs_len=ZSTD_compressBound(huffman_bs_len);
        zstd_bs[block_id].resize(max_zstd_bs_len);
        uint64_t zstd_bs_len=ZSTD_compress(zstd_bs[block_id].data(),max_zstd_bs_len,huffman_bs_vec.data(),huffman_bs_len,3);
        zstd_bs[block_id].resize(zstd_bs_len);
    }
    double end_time=omp_get_wtime();
    delete[] corner_qb;
    delete[] err_qb;
    delete[] no_corner_err_qb;
    printf("Compression time= %f seconds\n",end_time-start_time);
    printf("Bandwidth per cpu= %f MB/s\n",params.element_num_to_process*header.element_size*header.word_size/((end_time-start_time)*1024*1024*available_cpus));
    uint64_t total_compressed_size=0;
    for (auto block_id:params.blocks_to_process)
        total_compressed_size+=zstd_bs[block_id].size();
    printf("Total compressed size= %lu bytes\n",total_compressed_size);
    double br=static_cast<double>(total_compressed_size*8)/(params.element_num_to_process*header.element_size);
    double cr=static_cast<double>(header.word_size*8)/br;
    printf("BR= %f bits/value, CR= %f\n",br,cr);
}