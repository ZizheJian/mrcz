#pragma once

namespace mrcz
{
    template<typename T>
    void param_class<T>::decompress_irregular(string field)
    {
        uint32_t n0=header.element_shape[0];
        uint32_t n1=header.element_shape[1];
        uint32_t n2=header.element_shape[2];
        uint32_t n12=n1*n2;
        uint32_t n012=n0*n1*n2;
        vector<T> L0=L[0];
        vector<T> L1=L[1];
        vector<T> L2=L[2];
        T abs_eb=this->abs_eb[field];
        T coeff_eb=abs_eb*coeff_eb_scale;
        int64_t *coeff_qb=new int64_t[element_num_to_process*n012];
        int64_t *err_qb=new int64_t[element_num_to_process*n012];
        decompressed_data.resize(element_num_to_process*n012);
        double start_time=omp_get_wtime();
        #pragma omp parallel for schedule(static) num_threads(max_cpus)
        for (uint64_t block_id=0;block_id<blocks.size();block_id++)
        {
            uint64_t begin_index=blocks[block_id].begin_index;
            uint64_t end_index=blocks[block_id].end_index;
            // 无损解压
            uint64_t huffman_bs_len=ZSTD_getFrameContentSize(blocks[block_id].zstd_bs.data(),blocks[block_id].zstd_bs.size());
            vector<unsigned char> huffman_bs_vec(huffman_bs_len);
            if (ZSTD_decompress(huffman_bs_vec.data(),huffman_bs_len,blocks[block_id].zstd_bs.data(),blocks[block_id].zstd_bs.size())!=huffman_bs_len)
                throw runtime_error("ZSTD_decompress error in decompress()");
            const unsigned char *huffman_bs=huffman_bs_vec.data();
            if (if_zero_map)
                lossless_decompress_coeff_qb_with_zero_map(begin_index,end_index,coeff_qb,huffman_bs,n012);
            else
                lossless_decompress_coeff_qb(begin_index,end_index,coeff_qb,huffman_bs,n012);
            huffman_encoder_class<int64_t> err_encoder;
            err_encoder.load(huffman_bs);
            err_encoder.decode(huffman_bs,err_qb+n012*begin_index,n012*(end_index-begin_index));
            
            // 然后重构剩余的数据
            vector<T> temp_element1(n012,0),temp_element2(n012,0);
            for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
            {
                T *element=decompressed_data.data()+element_id*n012;
                int64_t *element_coeff_qb=coeff_qb+element_id*n012;
                int64_t *element_err_qb=err_qb+element_id*n012;
                // 重构coeff
                for (uint32_t i=0;i<n012;i++)
                    temp_element1[i]=element_coeff_qb[i]*2*coeff_eb;
                // 反变换
                inverse_transform_irregular(temp_element1,temp_element2,n0,n1,n2,n12,n012,L0,L1,L2);
                for (uint32_t j=0;j<n012;j++)
                    element[j]+=(temp_element2[j]-element_err_qb[j]*2*abs_eb);
            }
        }
        double end_time=omp_get_wtime();
        delete[] err_qb;
        delete[] coeff_qb;
        printf("Decompression time= %f seconds\n",end_time-start_time);
        printf("\033[91mBandwidth per cpu= %f MB/s\033[0m\n",element_num_to_process*n012*header.word_size/((end_time-start_time)*1024*1024*max_cpus));
    }
}