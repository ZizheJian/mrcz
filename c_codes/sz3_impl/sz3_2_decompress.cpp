#pragma once

#include "../mrcz_core/mrcz_param.hpp"

template<typename T>
void sz3_2_decompress(mrcz::param_class<T> &params,string field)
{
	auto &header=params.header;
    auto &blocks=params.blocks;
    auto &decompressed_data=params.decompressed_data;
    auto &x=params.x;

    uint32_t n=header.element_shape[0];
    uint32_t n2=n*n;
    uint32_t n3=n2*n;
    T abs_eb=params.abs_eb[field];
    int64_t *corner_qb=new int64_t[8*params.element_num_to_process];
    int64_t *err_qb=new int64_t[params.element_num_to_process*header.element_size];
    int64_t *no_corner_err_qb=new int64_t[params.element_num_to_process*(header.element_size-8)];
    vector<vector<uint32_t>> left_indices,right_indices,mid_indices,compressed_indices;
    vector<vector<T>> left_weight,right_weight;
    decompressed_data.resize(params.element_num_to_process*header.element_size);
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
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
    for (uint64_t block_id=0;block_id<blocks.size();block_id++) //这个block_id是假的，指的是在需要处理的blocks中的编号，而不是header文件中的编号
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
        err_encoder.decode(huffman_bs,no_corner_err_qb+(header.element_size-8)*begin_index,(header.element_size-8)*(end_index-begin_index));
        // 重构每个element角上的数据
        for (uint64_t i=begin_index;i<end_index;i++)
        {
            T *element=decompressed_data.data()+i*header.element_size;
            int64_t *element_corner_qb=corner_qb+8*i;
            element[0]=element_corner_qb[0]*2*abs_eb;
            element[n-1]=element_corner_qb[1]*2*abs_eb;
            element[(n-1)*n]=element_corner_qb[2]*2*abs_eb;
            element[n*n-1]=element_corner_qb[3]*2*abs_eb;
            element[(n-1)*n*n]=element_corner_qb[4]*2*abs_eb;
            element[(n-1)*(n*n+1)]=element_corner_qb[5]*2*abs_eb;
            element[(n-1)*n*(n+1)]=element_corner_qb[6]*2*abs_eb;
            element[n*n*n-1]=element_corner_qb[7]*2*abs_eb;
        }
        //把no_corner_err_qb转换成err_qb
        for (uint64_t i=begin_index;i<end_index;i++)
        {
            int64_t *element_err_qb=err_qb+i*header.element_size;
            int64_t *element_no_corner_err_qb=no_corner_err_qb+i*(header.element_size-8);
            memcpy(element_err_qb+1,element_no_corner_err_qb,sizeof(int64_t)*(n-2));
            memcpy(element_err_qb+n,element_no_corner_err_qb+n-2,sizeof(int64_t)*(n2-2*n));
            memcpy(element_err_qb+n2-n+1,element_no_corner_err_qb+n2-n-2,sizeof(int64_t)*(n-2));
            memcpy(element_err_qb+n2,element_no_corner_err_qb+n2-4,sizeof(int64_t)*(n3-2*n2));
            memcpy(element_err_qb+n3-n2+1,element_no_corner_err_qb+n3-n2-4,sizeof(int64_t)*(n-2));
            memcpy(element_err_qb+n3-n2+n,element_no_corner_err_qb+n3-n2+n-6,sizeof(int64_t)*(n2-2*n));
            memcpy(element_err_qb+n3-n+1,element_no_corner_err_qb+n3-n-6,sizeof(int64_t)*(n-2));
        }
        // 然后重构剩余的数据
        for (uint64_t i=begin_index;i<end_index;i++)
        {
            T *element=decompressed_data.data()+i*header.element_size;
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
                            element[mid_pos]=pred_value-2*abs_eb*element_err_qb[mid_pos];
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
                            element[mid_pos]=pred_value-2*abs_eb*element_err_qb[mid_pos];
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
                            element[mid_pos]=pred_value-2*abs_eb*element_err_qb[mid_pos];
                        }
            }
        }
    }
    double end_time=omp_get_wtime();
    delete[] corner_qb;
    delete[] err_qb;
    delete[] no_corner_err_qb;
    printf("Decompression time= %f seconds\n",end_time-start_time);
    printf("\033[91mBandwidth per cpu= %f MB/s\n\033[0m",params.element_num_to_process*header.element_size*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
}