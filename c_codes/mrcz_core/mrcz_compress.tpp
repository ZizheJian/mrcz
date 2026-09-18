#pragma once

namespace mrcz
{
    template<typename T>
    template<bool measure_phase_time,bool ablation>
    void param_class<T>::compress_mrcz(string field)
    {
        uint32_t n0=header.c_element_shape[0];
        uint32_t n1=header.c_element_shape[1];
        uint32_t n2=header.c_element_shape[2];
        uint32_t n12=n1*n2;
        uint32_t n012=n0*n1*n2;
        auto &L0=L[0];
        auto &L1=L[1];
        auto &L2=L[2];
        auto &neg_L0=neg_L[0];
        auto &neg_L1=neg_L[1];
        auto &neg_L2=neg_L[2];
        T abs_eb=this->abs_eb[field];
        T coeff_eb=abs_eb*coeff_eb_scale;
        T neg_2_coeff_eb=1/(2*coeff_eb);
        T inv_2_eb=1/(2*abs_eb);
        vector<double> forward_time(measure_phase_time?partitions.size():0,0),coeff_quant_time(measure_phase_time?partitions.size():0,0),inverse_time(measure_phase_time?partitions.size():0,0),residual_quant_time(measure_phase_time?partitions.size():0,0),coeff_huffman_time(measure_phase_time?partitions.size():0,0),residual_huffman_time(measure_phase_time?partitions.size():0,0),zstd_time(measure_phase_time?partitions.size():0,0);
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
            auto &block=partitions[block_id];
            uint64_t begin_index=block.begin_index;
            uint64_t end_index=block.end_index;
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
            vector<T> temp_element1(n012,0),temp_element2(n012,0);
            for (uint64_t element_id=begin_index;element_id<end_index;element_id++)
            {
                T *element=data.data()+element_id*n012;
                uint64_t local_element_id=element_id-begin_index;
                int64_t *element_coeff_qb=coeff_qb+local_element_id*n012;
                int64_t *element_err_qb=err_qb+local_element_id*n012;
                if constexpr(measure_phase_time) phase_start_time=omp_get_wtime();
                forward_transform(temp_element1,temp_element2,element,n0,n1,n2,n12,n012,L0,L1,L2);
                if constexpr(measure_phase_time) forward_time[block_id]+=omp_get_wtime()-phase_start_time;
                if constexpr(measure_phase_time) phase_start_time=omp_get_wtime();
                quantize_and_overwrite_array(temp_element1.data(),element_coeff_qb,n012,coeff_eb,neg_2_coeff_eb);
                if constexpr(measure_phase_time) coeff_quant_time[block_id]+=omp_get_wtime()-phase_start_time;
                if constexpr(measure_phase_time) phase_start_time=omp_get_wtime();
                inverse_transform(temp_element1,temp_element2,n0,n1,n2,n12,n012,neg_L0,neg_L1,neg_L2);
                if constexpr(measure_phase_time) inverse_time[block_id]+=omp_get_wtime()-phase_start_time;
                if constexpr(measure_phase_time) phase_start_time=omp_get_wtime();
                quantize_array(temp_element2.data(),element,element_err_qb,n012,abs_eb,inv_2_eb);
                if constexpr(measure_phase_time) residual_quant_time[block_id]+=omp_get_wtime()-phase_start_time;
            }
            // 无损压缩
            vector<unsigned char> huffman_bs_vec(2*n012*(end_index-begin_index)*sizeof(int64_t),0);
            unsigned char *huffman_bs=huffman_bs_vec.data();
            if constexpr(measure_phase_time) phase_start_time=omp_get_wtime();
            lossless_compress_coeff_qb<ablation>(0,block_element_num,coeff_qb,huffman_bs,n012);
            if constexpr(measure_phase_time) coeff_huffman_time[block_id]=omp_get_wtime()-phase_start_time;
            if constexpr(measure_phase_time) phase_start_time=omp_get_wtime();
            huffman_encoder_class<int64_t> err_encoder;
            err_encoder.preprocess_encode(err_qb,block_value_num);
            err_encoder.save(huffman_bs);
            err_encoder.encode(err_qb,block_value_num,huffman_bs);
            if constexpr(measure_phase_time) residual_huffman_time[block_id]=omp_get_wtime()-phase_start_time;
            uint64_t huffman_bs_len=huffman_bs-huffman_bs_vec.data();
            uint64_t max_zstd_bs_len=ZSTD_compressBound(huffman_bs_len);
            block.zstd_bs.resize(max_zstd_bs_len);
            if constexpr(measure_phase_time) phase_start_time=omp_get_wtime();
            uint64_t zstd_bs_len=ZSTD_compress(block.zstd_bs.data(),max_zstd_bs_len,huffman_bs_vec.data(),huffman_bs_len,3);
            if constexpr(measure_phase_time) zstd_time[block_id]=omp_get_wtime()-phase_start_time;
            block.zstd_bs.resize(zstd_bs_len);
        }
        for (int thread_id=0;thread_id<max_cpus;thread_id++)
        {
            delete[] coeff_qb_by_thread[thread_id];
            delete[] err_qb_by_thread[thread_id];
        }
        double end_time=omp_get_wtime();
        if constexpr(measure_phase_time)
        {
            double total_forward_time=0,total_coeff_quant_time=0,total_inverse_time=0,total_residual_quant_time=0,total_coeff_huffman_time=0,total_residual_huffman_time=0,total_zstd_time=0;
            for (uint64_t block_id=0;block_id<partitions.size();block_id++)
            {
                total_forward_time+=forward_time[block_id];
                total_coeff_quant_time+=coeff_quant_time[block_id];
                total_inverse_time+=inverse_time[block_id];
                total_residual_quant_time+=residual_quant_time[block_id];
                total_coeff_huffman_time+=coeff_huffman_time[block_id];
                total_residual_huffman_time+=residual_huffman_time[block_id];
                total_zstd_time+=zstd_time[block_id];
            }
            double measured_time=total_forward_time+total_coeff_quant_time+total_inverse_time+total_residual_quant_time+total_coeff_huffman_time+total_residual_huffman_time+total_zstd_time;
            printf("MRCZ compression accumulated partition CPU times: forward= %f, coefficient quantization= %f, inverse= %f, residual quantization= %f, coefficient Huffman= %f, residual Huffman= %f, Zstd= %f seconds\n",total_forward_time,total_coeff_quant_time,total_inverse_time,total_residual_quant_time,total_coeff_huffman_time,total_residual_huffman_time,total_zstd_time);
            if (measured_time>0)
                printf("MRCZ compression stage shares: forward= %.2f%%, coefficient quantization= %.2f%%, inverse= %.2f%%, residual quantization= %.2f%%, coefficient Huffman= %.2f%%, residual Huffman= %.2f%%, Zstd= %.2f%%\n",100*total_forward_time/measured_time,100*total_coeff_quant_time/measured_time,100*total_inverse_time/measured_time,100*total_residual_quant_time/measured_time,100*total_coeff_huffman_time/measured_time,100*total_residual_huffman_time/measured_time,100*total_zstd_time/measured_time);
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
