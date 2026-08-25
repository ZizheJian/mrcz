#pragma once

namespace mrcz
{
    template<typename T>
    void param_class<T>::compress_mrcz(string field)
    {
        uint32_t n0=header.element_shape[0];
        uint32_t n1=header.element_shape[1];
        uint32_t n2=header.element_shape[2];
        uint32_t n12=n1*n2;
        uint32_t n012=n0*n1*n2;
        vector<T> L0=L[0];
        vector<T> L1=L[1];
        vector<T> L2=L[2];
        vector<T> neg_L0=neg_L[0];
        vector<T> neg_L1=neg_L[1];
        vector<T> neg_L2=neg_L[2];
        T abs_eb=this->abs_eb[field];
        T coeff_eb=abs_eb*coeff_eb_scale;
        T neg_2_coeff_eb=1/(2*coeff_eb);
        T neg_2_eb=1/(2*abs_eb);
        bool element_batched=method=="MRCZ-3D-GLL-element-batched" or method=="MRCZ-3D-GLL-element-batched-original-transpose";
        bool direct_coefficient_major=method=="MRCZ-3D-GLL-element-batched";
        uint32_t batch_width=is_same_v<T,float>?8:4;
        int64_t *coeff_qb=new int64_t[element_num_to_process*n012];
        int64_t *err_qb=new int64_t[element_num_to_process*n012];
        char *profile_env=getenv("MRCZ_PROFILE_PHASES");
        bool profile_phases=profile_env!=nullptr and string(profile_env)!="0";
        uint64_t profile_block_num=profile_phases?partitions.size():0;
        vector<double> packing_time_by_block(profile_block_num,0),forward_time_by_block(profile_block_num,0),coeff_quant_time_by_block(profile_block_num,0),inverse_time_by_block(profile_block_num,0),residual_quant_time_by_block(profile_block_num,0);
        vector<double> coeff_huffman_time_by_block(profile_block_num,0),residual_huffman_time_by_block(profile_block_num,0),zstd_time_by_block(profile_block_num,0);
        double start_time=omp_get_wtime();
        #pragma omp parallel for schedule(static) num_threads(max_cpus)
        for (uint64_t block_id=0;block_id<partitions.size();block_id++)
        {
            auto &block=partitions[block_id];
            uint64_t begin_index=block.begin_index;
            uint64_t end_index=block.end_index;
            uint64_t block_element_num=end_index-begin_index;
            double phase_start_time=0;
            // 然后量化剩余的数据
            if (element_batched)
            {
                vector<T> element_batch(n012*batch_width,0),temp_element1(n012*batch_width,0),temp_element2(n012*batch_width,0);
                int64_t *block_coeff_qb=coeff_qb+begin_index*n012;
                uint64_t full_batch_end=begin_index+(block_element_num/batch_width)*batch_width;
                for (uint64_t batch_begin=begin_index;batch_begin<full_batch_end;batch_begin+=batch_width)
                {
                    if (profile_phases) phase_start_time=omp_get_wtime();
                    for (uint32_t lane=0;lane<batch_width;lane++)
                    {
                        T *element=data.data()+(batch_begin+lane)*n012;
                        for (uint32_t pos=0;pos<n012;pos++)
                            element_batch[pos*batch_width+lane]=element[pos];
                    }
                    if (profile_phases) packing_time_by_block[block_id]+=omp_get_wtime()-phase_start_time;
                    if (profile_phases) phase_start_time=omp_get_wtime();
                    forward_transform_element_batched(temp_element1,temp_element2,element_batch.data(),batch_width,n0,n1,n2,n12,n012,L0,L1,L2);
                    if (profile_phases) forward_time_by_block[block_id]+=omp_get_wtime()-phase_start_time;
                    if (profile_phases) phase_start_time=omp_get_wtime();
                    uint64_t local_batch_begin=batch_begin-begin_index;
                    for (uint32_t pos=0;pos<n012;pos++)
                        for (uint32_t lane=0;lane<batch_width;lane++)
                        {
                            uint64_t batch_pos=pos*batch_width+lane;
                            int64_t qb=static_cast<int64_t>(round(temp_element1[batch_pos]*neg_2_coeff_eb));
                            if (direct_coefficient_major)
                                block_coeff_qb[pos*block_element_num+local_batch_begin+lane]=qb;
                            else
                                coeff_qb[(batch_begin+lane)*n012+pos]=qb;
                            temp_element1[batch_pos]=qb*2*coeff_eb;
                        }
                    if (profile_phases) coeff_quant_time_by_block[block_id]+=omp_get_wtime()-phase_start_time;
                    if (profile_phases) phase_start_time=omp_get_wtime();
                    inverse_transform_element_batched(temp_element1,temp_element2,batch_width,n0,n1,n2,n12,n012,neg_L0,neg_L1,neg_L2);
                    if (profile_phases) inverse_time_by_block[block_id]+=omp_get_wtime()-phase_start_time;
                    if (profile_phases) phase_start_time=omp_get_wtime();
                    for (uint32_t lane=0;lane<batch_width;lane++)
                    {
                        T *element=data.data()+(batch_begin+lane)*n012;
                        int64_t *element_err_qb=err_qb+(batch_begin+lane)*n012;
                        for (uint32_t pos=0;pos<n012;pos++)
                            element_err_qb[pos]=static_cast<int64_t>(round((temp_element2[pos*batch_width+lane]-element[pos])*neg_2_eb));
                    }
                    if (profile_phases) residual_quant_time_by_block[block_id]+=omp_get_wtime()-phase_start_time;
                }
                vector<T> scalar_temp_element1(n012,0),scalar_temp_element2(n012,0);
                for (uint64_t element_id=full_batch_end;element_id<end_index;element_id++)
                {
                    T *element=data.data()+element_id*n012;
                    int64_t *element_err_qb=err_qb+element_id*n012;
                    if (profile_phases) phase_start_time=omp_get_wtime();
                    forward_transform(scalar_temp_element1,scalar_temp_element2,element,n0,n1,n2,n12,n012,L0,L1,L2);
                    if (profile_phases) forward_time_by_block[block_id]+=omp_get_wtime()-phase_start_time;
                    if (profile_phases) phase_start_time=omp_get_wtime();
                    uint64_t local_element_id=element_id-begin_index;
                    for (uint32_t pos=0;pos<n012;pos++)
                    {
                        int64_t qb=static_cast<int64_t>(round(scalar_temp_element1[pos]*neg_2_coeff_eb));
                        if (direct_coefficient_major)
                            block_coeff_qb[pos*block_element_num+local_element_id]=qb;
                        else
                            coeff_qb[element_id*n012+pos]=qb;
                        scalar_temp_element1[pos]=qb*2*coeff_eb;
                    }
                    if (profile_phases) coeff_quant_time_by_block[block_id]+=omp_get_wtime()-phase_start_time;
                    if (profile_phases) phase_start_time=omp_get_wtime();
                    inverse_transform(scalar_temp_element1,scalar_temp_element2,n0,n1,n2,n12,n012,neg_L0,neg_L1,neg_L2);
                    if (profile_phases) inverse_time_by_block[block_id]+=omp_get_wtime()-phase_start_time;
                    if (profile_phases) phase_start_time=omp_get_wtime();
                    quantize_err_array(scalar_temp_element2.data(),element,element_err_qb,n012,abs_eb,neg_2_eb);
                    if (profile_phases) residual_quant_time_by_block[block_id]+=omp_get_wtime()-phase_start_time;
                }
            }
            else
            {
                vector<T> temp_element1(n012,0),temp_element2(n012,0);
                for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
                {
                    T *element=data.data()+element_id*n012;
                    int64_t *element_f_qb=coeff_qb+element_id*n012;
                    int64_t *element_err_qb=err_qb+element_id*n012;
                    if (profile_phases) phase_start_time=omp_get_wtime();
                    forward_transform(temp_element1,temp_element2,element,n0,n1,n2,n12,n012,L0,L1,L2);
                    if (profile_phases) forward_time_by_block[block_id]+=omp_get_wtime()-phase_start_time;
                    if (profile_phases) phase_start_time=omp_get_wtime();
                    quantize_and_overwrite_array1(temp_element1.data(),element_f_qb,n012,coeff_eb,neg_2_coeff_eb);
                    if (profile_phases) coeff_quant_time_by_block[block_id]+=omp_get_wtime()-phase_start_time;
                    if (profile_phases) phase_start_time=omp_get_wtime();
                    inverse_transform(temp_element1,temp_element2,n0,n1,n2,n12,n012,neg_L0,neg_L1,neg_L2);
                    if (profile_phases) inverse_time_by_block[block_id]+=omp_get_wtime()-phase_start_time;
                    if (profile_phases) phase_start_time=omp_get_wtime();
                    quantize_err_array(temp_element2.data(),element,element_err_qb,n012,abs_eb,neg_2_eb);
                    if (profile_phases) residual_quant_time_by_block[block_id]+=omp_get_wtime()-phase_start_time;
                }
            }
            // 无损压缩
            vector<unsigned char> huffman_bs_vec(2*n012*(end_index-begin_index)*sizeof(int64_t),0);
            unsigned char *huffman_bs=huffman_bs_vec.data();
            if (profile_phases) phase_start_time=omp_get_wtime();
            if (direct_coefficient_major)
            {
                huffman_encoder_class<int64_t> coeff_encoder;
                coeff_encoder.preprocess_encode(coeff_qb+begin_index*n012,n012*block_element_num);
                coeff_encoder.save(huffman_bs);
                coeff_encoder.encode(coeff_qb+begin_index*n012,n012*block_element_num,huffman_bs);
            }
            else
                lossless_compress_coeff_qb(begin_index,end_index,coeff_qb,huffman_bs,n012);
            if (profile_phases) coeff_huffman_time_by_block[block_id]=omp_get_wtime()-phase_start_time;
            if (profile_phases) phase_start_time=omp_get_wtime();
            huffman_encoder_class<int64_t> err_encoder;
            err_encoder.preprocess_encode(err_qb+begin_index*n012,n012*(end_index-begin_index));
            err_encoder.save(huffman_bs);
            err_encoder.encode(err_qb+begin_index*n012,n012*(end_index-begin_index),huffman_bs);
            if (profile_phases) residual_huffman_time_by_block[block_id]=omp_get_wtime()-phase_start_time;
            uint64_t huffman_bs_len=huffman_bs-huffman_bs_vec.data();
            uint64_t max_zstd_bs_len=ZSTD_compressBound(huffman_bs_len);
            block.zstd_bs.resize(max_zstd_bs_len);
            if (profile_phases) phase_start_time=omp_get_wtime();
            uint64_t zstd_bs_len=ZSTD_compress(block.zstd_bs.data(),max_zstd_bs_len,huffman_bs_vec.data(),huffman_bs_len,3);
            if (profile_phases) zstd_time_by_block[block_id]=omp_get_wtime()-phase_start_time;
            block.zstd_bs.resize(zstd_bs_len);

        }
        double end_time=omp_get_wtime();
        delete[] coeff_qb;
        delete[] err_qb;
        if (profile_phases)
        {
            double packing_time=0,forward_time=0,coeff_quant_time=0,inverse_time=0,residual_quant_time=0,coeff_huffman_time=0,residual_huffman_time=0,zstd_time=0;
            for (uint64_t block_id=0;block_id<profile_block_num;block_id++)
            {
                packing_time+=packing_time_by_block[block_id];
                forward_time+=forward_time_by_block[block_id];
                coeff_quant_time+=coeff_quant_time_by_block[block_id];
                inverse_time+=inverse_time_by_block[block_id];
                residual_quant_time+=residual_quant_time_by_block[block_id];
                coeff_huffman_time+=coeff_huffman_time_by_block[block_id];
                residual_huffman_time+=residual_huffman_time_by_block[block_id];
                zstd_time+=zstd_time_by_block[block_id];
            }
            double measured_time=packing_time+forward_time+coeff_quant_time+inverse_time+residual_quant_time+coeff_huffman_time+residual_huffman_time+zstd_time;
            printf("MRCZ compression accumulated partition CPU times: element packing= %f, forward= %f, coefficient quantization= %f, inverse= %f, residual quantization= %f, coefficient Huffman= %f, residual Huffman= %f, Zstd= %f seconds\n",packing_time,forward_time,coeff_quant_time,inverse_time,residual_quant_time,coeff_huffman_time,residual_huffman_time,zstd_time);
            if (measured_time>0)
                printf("MRCZ compression stage shares: element packing= %.2f%%, forward= %.2f%%, coefficient quantization= %.2f%%, inverse= %.2f%%, residual quantization= %.2f%%, coefficient Huffman= %.2f%%, residual Huffman= %.2f%%, Zstd= %.2f%%\n",100*packing_time/measured_time,100*forward_time/measured_time,100*coeff_quant_time/measured_time,100*inverse_time/measured_time,100*residual_quant_time/measured_time,100*coeff_huffman_time/measured_time,100*residual_huffman_time/measured_time,100*zstd_time/measured_time);
        }
        printf("Compression time= %f seconds\n",end_time-start_time);
        printf("\033[91mField= %s, compression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),element_num_to_process*n012*header.word_size/((end_time-start_time)*1024*1024*max_cpus));
    }

    template<typename T>
    inline void param_class<T>::forward_transform(vector<T> &temp_element1,vector<T> &temp_element2,T *element,uint32_t n0,uint32_t n1,uint32_t n2,uint32_t n12,uint32_t n012,
        vector<T> &L0,vector<T> &L1,vector<T> &L2)
    {
        T *temp1=temp_element1.data(),*temp2=temp_element2.data();
        T *l0=L0.data(),*l1=L1.data(),*l2=L2.data();
        memset(temp1,0,sizeof(T)*n012);
        for (uint32_t a=0;a<n0;a++)
        {
            T *input=element+a*n12;
            for (uint32_t i=0;i<n0;i++)
            {
                T *output=temp1+i*n12;
                T weight=l0[a*n0+i];
                for (uint32_t bn2pc=0;bn2pc<n1*n2;bn2pc++)
                    output[bn2pc]+=input[bn2pc]*weight;
            }
        }
        memset(temp2,0,sizeof(T)*n012);
        for (uint32_t in12=0;in12<n012;in12+=n12)
            for (uint32_t b=0;b<n1;b++)
            {
                T *input=temp1+in12+b*n2;
                for (uint32_t j=0;j<n1;j++)
                {
                    T *output=temp2+in12+j*n2;
                    T weight=l1[b*n1+j];
                    for (uint32_t c=0;c<n2;c++)
                        output[c]+=input[c]*weight;
                }
            }
        memset(temp1,0,sizeof(T)*n012);
        for (uint32_t in12pjn2=0;in12pjn2<n012;in12pjn2+=n2)
            for (uint32_t c=0;c<n2;c++)
            {
                T input=temp2[in12pjn2+c];
                T *output=temp1+in12pjn2;
                T *weights=l2+c*n2;
                for (uint32_t k=0;k<n2;k++)
                    output[k]+=input*weights[k];
            }
    }

    template<typename T>
    inline void param_class<T>::inverse_transform(vector<T> &temp_element1,vector<T> &temp_element2,uint32_t n0,uint32_t n1,uint32_t n2,uint32_t n12,uint32_t n012,
        vector<T> &neg_L0,vector<T> &neg_L1,vector<T> &neg_L2)
    {
        T *temp1=temp_element1.data(),*temp2=temp_element2.data();
        T *l0=neg_L0.data(),*l1=neg_L1.data(),*l2=neg_L2.data();
        memset(temp2,0,sizeof(T)*n012);
        for (uint32_t in12pjn2=0;in12pjn2<n012;in12pjn2+=n2)
            for (uint32_t k=0;k<n2;k++)
            {
                T input=temp1[in12pjn2+k];
                T *output=temp2+in12pjn2;
                T *weights=l2+k*n2;
                for (uint32_t c=0;c<n2;c++)
                    output[c]+=input*weights[c];
            }
        memset(temp1,0,sizeof(T)*n012);
        for (uint32_t in12=0;in12<n012;in12+=n12)
            for (uint32_t j=0;j<n1;j++)
            {
                T *input=temp2+in12+j*n2;
                for (uint32_t b=0;b<n1;b++)
                {
                    T *output=temp1+in12+b*n2;
                    T weight=l1[j*n1+b];
                    for (uint32_t c=0;c<n2;c++)
                        output[c]+=input[c]*weight;
                }
            }
        memset(temp2,0,sizeof(T)*n012);
        for (uint32_t i=0;i<n0;i++)
        {
            T *input=temp1+i*n12;
            for (uint32_t a=0;a<n0;a++)
            {
                T *output=temp2+a*n12;
                T weight=l0[i*n0+a];
                for (uint32_t bn2pc=0;bn2pc<n12;bn2pc++)
                    output[bn2pc]+=input[bn2pc]*weight;
            }
        }
    }
}
