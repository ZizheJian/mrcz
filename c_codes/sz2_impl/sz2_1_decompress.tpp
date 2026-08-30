#pragma once


#include "../mrcz_core/mrcz.hpp"
#include "../mrcz_core/huffman.hpp"
#include "../third_party/zstd/lib/zstd.h"

template<typename T>
void sz2_1_decompress(mrcz::param_class<T> &params,string field)
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
    uint32_t n012=header.element_size;
    T abs_eb=params.abs_eb[field];
    T coeff_eb=0.1*abs_eb;
    int64_t *err_qb=new int64_t[params.element_num_to_process*(n012+4)];
    decompressed_data.resize(params.element_num_to_process*n012);
    double start_time=omp_get_wtime();
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
    for (uint64_t block_id=0;block_id<blocks.size();block_id++) //这个block_id是假的，指的是在需要处理的blocks中的编号，而不是header文件中的编号
    {
        auto &block=blocks[block_id];
        uint64_t begin_index=block.begin_index;
        uint64_t end_index=block.end_index;
        int64_t *a_address=err_qb+begin_index*(n012+4);
        int64_t *b_address=a_address+end_index-begin_index;
        int64_t *c_address=b_address+end_index-begin_index;
        int64_t *d_address=c_address+end_index-begin_index;
        int64_t *block_err_qb=d_address+end_index-begin_index;
        // 先无损解压
        uint64_t huffman_bs_len=ZSTD_getFrameContentSize(block.zstd_bs.data(),block.zstd_bs.size());
        vector<unsigned char> huffman_bs_vec(huffman_bs_len);
        if (ZSTD_decompress(huffman_bs_vec.data(),huffman_bs_len,block.zstd_bs.data(),block.zstd_bs.size())!=huffman_bs_len)
            throw runtime_error("ZSTD_decompress error in decompress()");
        const unsigned char *huffman_bs=huffman_bs_vec.data();
        mrcz::huffman_encoder_class<int64_t> err_encoder;
        err_encoder.load(huffman_bs);
        err_encoder.decode(huffman_bs,err_qb+(n012+4)*begin_index,(n012+4)*(end_index-begin_index));

        for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
        {
            T *element=decompressed_data.data()+element_id*n012;
            T a=a_address[element_id-begin_index]*coeff_eb*2;
            T b=b_address[element_id-begin_index]*coeff_eb*2;
            T c=c_address[element_id-begin_index]*coeff_eb*2;
            T d=d_address[element_id-begin_index]*coeff_eb*2;
            vector<T> new_element(n012,0);
            for (uint32_t i0=0;i0<n0;i0++)
                for (uint32_t i1=0;i1<n1;i1++)
                    for (uint32_t i2=0;i2<n2;i2++)
                        element[i0*n12+i1*n2+i2]=a*x0[i0]+b*x1[i1]+c*x2[i2]+d-block_err_qb[(element_id-begin_index)*n012+i0*n12+i1*n2+i2]*2*abs_eb;
        }
    }
    double end_time=omp_get_wtime();
    delete[] err_qb;
    printf("Decompression time= %f seconds\n",end_time-start_time);
    printf("\033[91mField= %s, decompression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),params.element_num_to_process*n012*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
}
