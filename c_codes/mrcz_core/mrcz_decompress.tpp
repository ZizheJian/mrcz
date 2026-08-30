#pragma once

namespace mrcz
{
    template<typename T>
    void param_class<T>::decompress_mrcz(string field)
    {
        uint32_t n0=header.c_element_shape[0];
        uint32_t n1=header.c_element_shape[1];
        uint32_t n2=header.c_element_shape[2];
        uint32_t n12=n1*n2;
        uint32_t n012=n0*n1*n2;
        vector<T> neg_L0=neg_L[0];
        vector<T> neg_L1=neg_L[1];
        vector<T> neg_L2=neg_L[2];
        T abs_eb=this->abs_eb[field];
        T coeff_eb=abs_eb*coeff_eb_scale;
        bool element_batched=method=="MRCZ-element-batched" or method=="MRCZ-element-batched-original-transpose";
        bool direct_coefficient_major=method=="MRCZ-element-batched";
        uint32_t batch_width=is_same_v<T,float>?8:4;
        decompressed_data.resize(element_num_to_process*n012);
        char *profile_env=getenv("MRCZ_PROFILE_PHASES");
        bool profile_phases=profile_env!=nullptr and string(profile_env)!="0";
        uint64_t profile_block_num=profile_phases?partitions.size():0;
        vector<double> zstd_time_by_block(profile_block_num,0),coeff_huffman_time_by_block(profile_block_num,0),residual_huffman_time_by_block(profile_block_num,0);
        vector<double> coeff_reconstruct_time_by_block(profile_block_num,0),inverse_time_by_block(profile_block_num,0),residual_reconstruct_time_by_block(profile_block_num,0);
        vector<int64_t*> coeff_qb_by_thread(max_cpus,nullptr),err_qb_by_thread(max_cpus,nullptr);
        vector<uint64_t> qb_capacity_by_thread(max_cpus,0);
        double start_time=omp_get_wtime();
        #pragma omp parallel for schedule(static) num_threads(max_cpus)
        for (uint64_t block_id=0;block_id<partitions.size();block_id++)
        {
            int thread_id=omp_get_thread_num();
            int64_t *&coeff_qb=coeff_qb_by_thread[thread_id];
            int64_t *&err_qb=err_qb_by_thread[thread_id];
            uint64_t &qb_capacity=qb_capacity_by_thread[thread_id];
            uint64_t begin_index=partitions[block_id].begin_index;
            uint64_t end_index=partitions[block_id].end_index;
            uint64_t block_element_num=end_index-begin_index;
            uint64_t block_value_num=block_element_num*n012;
            if (block_value_num>qb_capacity)
            {
                delete[] coeff_qb;
                delete[] err_qb;
                coeff_qb=new int64_t[block_value_num];
                err_qb=new int64_t[block_value_num];
                qb_capacity=block_value_num;
            }
            double phase_start_time=0;
            // 无损解压
            uint64_t huffman_bs_len=ZSTD_getFrameContentSize(partitions[block_id].zstd_bs.data(),partitions[block_id].zstd_bs.size());
            vector<unsigned char> huffman_bs_vec(huffman_bs_len);
            if (profile_phases) phase_start_time=omp_get_wtime();
            if (ZSTD_decompress(huffman_bs_vec.data(),huffman_bs_len,partitions[block_id].zstd_bs.data(),partitions[block_id].zstd_bs.size())!=huffman_bs_len)
                throw runtime_error("ZSTD_decompress error in decompress()");
            if (profile_phases) zstd_time_by_block[block_id]=omp_get_wtime()-phase_start_time;
            const unsigned char *huffman_bs=huffman_bs_vec.data();
            if (profile_phases) phase_start_time=omp_get_wtime();
            if (direct_coefficient_major)
            {
                huffman_encoder_class<int64_t> coeff_encoder;
                coeff_encoder.load(huffman_bs);
                coeff_encoder.decode(huffman_bs,coeff_qb,block_value_num);
            }
            else
                lossless_decompress_coeff_qb(0,block_element_num,coeff_qb,huffman_bs,n012);
            if (profile_phases) coeff_huffman_time_by_block[block_id]=omp_get_wtime()-phase_start_time;
            if (profile_phases) phase_start_time=omp_get_wtime();
            huffman_encoder_class<int64_t> err_encoder;
            err_encoder.load(huffman_bs);
            err_encoder.decode(huffman_bs,err_qb,block_value_num);
            if (profile_phases) residual_huffman_time_by_block[block_id]=omp_get_wtime()-phase_start_time;
            
            // 然后重构剩余的数据
            if (element_batched)
            {
                vector<T> temp_element1(n012*batch_width,0),temp_element2(n012*batch_width,0);
                int64_t *block_coeff_qb=coeff_qb;
                uint64_t full_batch_end=begin_index+(block_element_num/batch_width)*batch_width;
                for (uint64_t batch_begin=begin_index;batch_begin<full_batch_end;batch_begin+=batch_width)
                {
                    if (profile_phases) phase_start_time=omp_get_wtime();
                    uint64_t local_batch_begin=batch_begin-begin_index;
                    for (uint32_t pos=0;pos<n012;pos++)
                        for (uint32_t lane=0;lane<batch_width;lane++)
                            if (direct_coefficient_major)
                                temp_element1[pos*batch_width+lane]=block_coeff_qb[pos*block_element_num+local_batch_begin+lane]*2*coeff_eb;
                            else
                                temp_element1[pos*batch_width+lane]=coeff_qb[(local_batch_begin+lane)*n012+pos]*2*coeff_eb;
                    if (profile_phases) coeff_reconstruct_time_by_block[block_id]+=omp_get_wtime()-phase_start_time;
                    if (profile_phases) phase_start_time=omp_get_wtime();
                    inverse_transform_element_batched(temp_element1,temp_element2,batch_width,n0,n1,n2,n12,n012,neg_L0,neg_L1,neg_L2);
                    if (profile_phases) inverse_time_by_block[block_id]+=omp_get_wtime()-phase_start_time;
                    if (profile_phases) phase_start_time=omp_get_wtime();
                    for (uint32_t lane=0;lane<batch_width;lane++)
                    {
                        T *element=decompressed_data.data()+(batch_begin+lane)*n012;
                        int64_t *element_err_qb=err_qb+(local_batch_begin+lane)*n012;
                        for (uint32_t pos=0;pos<n012;pos++)
                            element[pos]+=(temp_element2[pos*batch_width+lane]-element_err_qb[pos]*2*abs_eb);
                    }
                    if (profile_phases) residual_reconstruct_time_by_block[block_id]+=omp_get_wtime()-phase_start_time;
                }
                vector<T> scalar_temp_element1(n012,0),scalar_temp_element2(n012,0);
                for (uint64_t element_id=full_batch_end;element_id<end_index;element_id++)
                {
                    T *element=decompressed_data.data()+element_id*n012;
                    uint64_t local_element_id=element_id-begin_index;
                    int64_t *element_err_qb=err_qb+local_element_id*n012;
                    if (profile_phases) phase_start_time=omp_get_wtime();
                    for (uint32_t pos=0;pos<n012;pos++)
                        if (direct_coefficient_major)
                            scalar_temp_element1[pos]=block_coeff_qb[pos*block_element_num+local_element_id]*2*coeff_eb;
                        else
                            scalar_temp_element1[pos]=coeff_qb[local_element_id*n012+pos]*2*coeff_eb;
                    if (profile_phases) coeff_reconstruct_time_by_block[block_id]+=omp_get_wtime()-phase_start_time;
                    if (profile_phases) phase_start_time=omp_get_wtime();
                    inverse_transform(scalar_temp_element1,scalar_temp_element2,n0,n1,n2,n12,n012,neg_L0,neg_L1,neg_L2);
                    if (profile_phases) inverse_time_by_block[block_id]+=omp_get_wtime()-phase_start_time;
                    if (profile_phases) phase_start_time=omp_get_wtime();
                    for (uint32_t pos=0;pos<n012;pos++)
                        element[pos]+=(scalar_temp_element2[pos]-element_err_qb[pos]*2*abs_eb);
                    if (profile_phases) residual_reconstruct_time_by_block[block_id]+=omp_get_wtime()-phase_start_time;
                }
            }
            else
            {
                vector<T> temp_element1(n012,0),temp_element2(n012,0);
                for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
                {
                    T *element=decompressed_data.data()+element_id*n012;
                    uint64_t local_element_id=element_id-begin_index;
                    int64_t *element_coeff_qb=coeff_qb+local_element_id*n012;
                    int64_t *element_err_qb=err_qb+local_element_id*n012;
                    if (profile_phases) phase_start_time=omp_get_wtime();
                    for (uint32_t i=0;i<n012;i++)
                        temp_element1[i]=element_coeff_qb[i]*2*coeff_eb;
                    if (profile_phases) coeff_reconstruct_time_by_block[block_id]+=omp_get_wtime()-phase_start_time;
                    if (profile_phases) phase_start_time=omp_get_wtime();
                    inverse_transform(temp_element1,temp_element2,n0,n1,n2,n12,n012,neg_L0,neg_L1,neg_L2);
                    if (profile_phases) inverse_time_by_block[block_id]+=omp_get_wtime()-phase_start_time;
                    if (profile_phases) phase_start_time=omp_get_wtime();
                    for (uint32_t j=0;j<n012;j++)
                        element[j]+=(temp_element2[j]-element_err_qb[j]*2*abs_eb);
                    if (profile_phases) residual_reconstruct_time_by_block[block_id]+=omp_get_wtime()-phase_start_time;
                }
            }
        }
        for (int thread_id=0;thread_id<max_cpus;thread_id++)
        {
            delete[] coeff_qb_by_thread[thread_id];
            delete[] err_qb_by_thread[thread_id];
        }
        double end_time=omp_get_wtime();
        if (profile_phases)
        {
            double zstd_time=0,coeff_huffman_time=0,residual_huffman_time=0,coeff_reconstruct_time=0,inverse_time=0,residual_reconstruct_time=0;
            for (uint64_t block_id=0;block_id<profile_block_num;block_id++)
            {
                zstd_time+=zstd_time_by_block[block_id];
                coeff_huffman_time+=coeff_huffman_time_by_block[block_id];
                residual_huffman_time+=residual_huffman_time_by_block[block_id];
                coeff_reconstruct_time+=coeff_reconstruct_time_by_block[block_id];
                inverse_time+=inverse_time_by_block[block_id];
                residual_reconstruct_time+=residual_reconstruct_time_by_block[block_id];
            }
            double measured_time=zstd_time+coeff_huffman_time+residual_huffman_time+coeff_reconstruct_time+inverse_time+residual_reconstruct_time;
            printf("MRCZ decompression accumulated partition CPU times: Zstd= %f, coefficient Huffman= %f, residual Huffman= %f, coefficient reconstruction= %f, inverse= %f, residual reconstruction= %f seconds\n",zstd_time,coeff_huffman_time,residual_huffman_time,coeff_reconstruct_time,inverse_time,residual_reconstruct_time);
            if (measured_time>0)
                printf("MRCZ decompression stage shares: Zstd= %.2f%%, coefficient Huffman= %.2f%%, residual Huffman= %.2f%%, coefficient reconstruction= %.2f%%, inverse= %.2f%%, residual reconstruction= %.2f%%\n",100*zstd_time/measured_time,100*coeff_huffman_time/measured_time,100*residual_huffman_time/measured_time,100*coeff_reconstruct_time/measured_time,100*inverse_time/measured_time,100*residual_reconstruct_time/measured_time);
        }
        printf("Decompression time= %f seconds\n",end_time-start_time);
        printf("\033[91mField= %s, decompression bandwidth per cpu= %f MB/s\033[0m\n",field.c_str(),element_num_to_process*n012*header.word_size/((end_time-start_time)*1024*1024*max_cpus));
    }
}
