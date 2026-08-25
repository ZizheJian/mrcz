#pragma once

# XXX 该文件保存2026-08-25完成的MRCZ coefficient multi-tree实验实现，不参与正式编译。
# XXX 双树按总阶数阈值0到4分组；三树按常数、线性、其余分组；四树按常数、线性、二次、其余分组。
# XXX 三个小数据上三树相对单树BR几何平均下降8.54%，四树下降9.67%。
# XXX jicf/S00的三个大partition样本上三树和四树仅分别下降1.13%和1.14%，四树相对三树仅下降0.012%。
# XXX 收益依赖数据区域且没有稳定速度优势，因此正式MRCZ恢复coefficient-major单树编码。

namespace mrcz::trash
{
    inline uint32_t coefficient_group_id(uint32_t order,uint32_t tree_num,uint32_t split_order)
    {
        if (tree_num==2)
            return order<=split_order?0:1;
        return min(order,tree_num-1);
    }

    template<typename T>
    void lossless_compress_coeff_qb_multitree(param_class<T> &params,uint64_t begin_index,uint64_t end_index,int64_t *coeff_qb,unsigned char *&huffman_bs,uint32_t n000,uint32_t tree_num,uint32_t split_order=1)
    {
        uint32_t n1=params.header.element_shape[1],n2=params.header.element_shape[2],n12=n1*n2;
        uint64_t element_num=end_index-begin_index;
        vector<uint32_t> group_mode_num(tree_num,0);
        for (uint32_t pos=0;pos<n000;pos++)
        {
            uint32_t i=pos/n12,j=(pos%n12)/n2,k=pos%n2,group_id=coefficient_group_id(i+j+k,tree_num,split_order);
            group_mode_num[group_id]++;
        }
        vector<vector<int64_t>> grouped_coeff_qb(tree_num);
        vector<uint64_t> group_index(tree_num,0);
        for (uint32_t group_id=0;group_id<tree_num;group_id++)
            grouped_coeff_qb[group_id].resize(group_mode_num[group_id]*element_num);
        for (uint32_t pos=0;pos<n000;pos++)
        {
            uint32_t i=pos/n12,j=(pos%n12)/n2,k=pos%n2,group_id=coefficient_group_id(i+j+k,tree_num,split_order);
            for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
                grouped_coeff_qb[group_id][group_index[group_id]++]=coeff_qb[element_id*n000+pos];
        }
        for (uint32_t group_id=0;group_id<tree_num;group_id++)
        {
            if (grouped_coeff_qb[group_id].empty())
                continue;
            huffman_encoder_class<int64_t> encoder;
            encoder.preprocess_encode(grouped_coeff_qb[group_id].data(),grouped_coeff_qb[group_id].size());
            encoder.save(huffman_bs);
            encoder.encode(grouped_coeff_qb[group_id].data(),grouped_coeff_qb[group_id].size(),huffman_bs);
        }
    }

    template<typename T>
    void lossless_decompress_coeff_qb_multitree(param_class<T> &params,uint64_t begin_index,uint64_t end_index,int64_t *coeff_qb,unsigned char const *&huffman_bs,uint32_t n000,uint32_t tree_num,uint32_t split_order=1)
    {
        uint32_t n1=params.header.element_shape[1],n2=params.header.element_shape[2],n12=n1*n2;
        uint64_t element_num=end_index-begin_index;
        vector<uint32_t> group_mode_num(tree_num,0);
        for (uint32_t pos=0;pos<n000;pos++)
        {
            uint32_t i=pos/n12,j=(pos%n12)/n2,k=pos%n2,group_id=coefficient_group_id(i+j+k,tree_num,split_order);
            group_mode_num[group_id]++;
        }
        vector<vector<int64_t>> grouped_coeff_qb(tree_num);
        for (uint32_t group_id=0;group_id<tree_num;group_id++)
        {
            grouped_coeff_qb[group_id].resize(group_mode_num[group_id]*element_num);
            if (grouped_coeff_qb[group_id].empty())
                continue;
            huffman_encoder_class<int64_t> encoder;
            encoder.load(huffman_bs);
            encoder.decode(huffman_bs,grouped_coeff_qb[group_id].data(),grouped_coeff_qb[group_id].size());
        }
        vector<uint64_t> group_index(tree_num,0);
        for (uint32_t pos=0;pos<n000;pos++)
        {
            uint32_t i=pos/n12,j=(pos%n12)/n2,k=pos%n2,group_id=coefficient_group_id(i+j+k,tree_num,split_order);
            for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
                coeff_qb[element_id*n000+pos]=grouped_coeff_qb[group_id][group_index[group_id]++];
        }
    }
}
