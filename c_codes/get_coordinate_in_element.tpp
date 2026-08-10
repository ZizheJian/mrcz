#pragma once

#include "mrcz_core/mrcz_param.hpp"

template<typename T>
void get_coordinate_in_element(mrcz::param_class<T> &params)
{
    vector<float> coordinate_float(params.element_num_to_process*3*params.header.element_size);
    ifstream data_file(params.data_path,ios::binary);
    data_file.seekg(params.field_begin_address["X"],ios::beg);
    data_file.read(reinterpret_cast<char*>(coordinate_float.data()),params.element_num_to_process*3*params.header.element_size*params.header.word_size);
    data_file.close();

    vector<T> coordinate(coordinate_float.begin(),coordinate_float.end());

    uint32_t n=params.header.element_shape[0];
    uint32_t n2=n*n;
    vector<T> segment_length(n-1,0);
    vector<T> segment_length_2(n-1,0);
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
    for (uint64_t block_id=0;block_id<params.partitions.size();block_id++)
    {
        vector<T> block_segment_length(n-1,0);
        vector<T> block_segment_length_2(n-1,0);
        auto &block=params.partitions[block_id];
        uint64_t begin_element_index=block.begin_index;
        uint64_t end_element_index=block.end_index;
        for (uint64_t element_id=begin_element_index;element_id<end_element_index;element_id++)
        {
            T *X=&coordinate[element_id*3*params.header.element_size];
            T *Y=&coordinate[element_id*3*params.header.element_size+params.header.element_size];
            T *Z=&coordinate[element_id*3*params.header.element_size+2*params.header.element_size];
            for (uint32_t i0=0;i0<n;i0++)
                for (uint32_t i1=0;i1<n;i1++)
                {
                    vector<T> line_segment_length(n-1,0);
                    T line_total_length=0;
                    for (uint32_t i2=0;i2<n-1;i2++)
                    {
                        T x1=X[i0*n2+i1*n+i2];
                        T y1=Y[i0*n2+i1*n+i2];
                        T z1=Z[i0*n2+i1*n+i2];
                        T x2=X[i0*n2+i1*n+i2+1];
                        T y2=Y[i0*n2+i1*n+i2+1];
                        T z2=Z[i0*n2+i1*n+i2+1];
                        line_segment_length[i2]=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
                        line_total_length+=line_segment_length[i2];
                    }
                    for (uint32_t i2=0;i2<n-1;i2++)
                    {
                        line_segment_length[i2]*=((T)2/line_total_length);
                        block_segment_length[i2]+=line_segment_length[i2];
                        block_segment_length_2[i2]+=line_segment_length[i2]*line_segment_length[i2];
                    }
                }
            for (uint32_t i0=0;i0<n;i0++)
                for (uint32_t i2=0;i2<n;i2++)
                {
                    vector<T> line_segment_length(n-1,0);
                    T line_total_length=0;
                    for (uint32_t i1=0;i1<n-1;i1++)
                    {
                        T x1=X[i0*n2+i1*n+i2];
                        T y1=Y[i0*n2+i1*n+i2];
                        T z1=Z[i0*n2+i1*n+i2];
                        T x2=X[i0*n2+(i1+1)*n+i2];
                        T y2=Y[i0*n2+(i1+1)*n+i2];
                        T z2=Z[i0*n2+(i1+1)*n+i2];
                        line_segment_length[i1]=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
                        line_total_length+=line_segment_length[i1];
                    }
                    for (uint32_t i1=0;i1<n-1;i1++)
                    {
                        line_segment_length[i1]*=((T)2/line_total_length);
                        block_segment_length[i1]+=line_segment_length[i1];
                        block_segment_length_2[i1]+=line_segment_length[i1]*line_segment_length[i1];
                    }
                }
            for (uint32_t i1=0;i1<n;i1++)
                for (uint32_t i2=0;i2<n;i2++)
                {
                    vector<T> line_segment_length(n-1,0);
                    T line_total_length=0;
                    for (uint32_t i0=0;i0<n-1;i0++)
                    {
                        T x1=X[i0*n2+i1*n+i2];
                        T y1=Y[i0*n2+i1*n+i2];
                        T z1=Z[i0*n2+i1*n+i2];
                        T x2=X[(i0+1)*n2+i1*n+i2];
                        T y2=Y[(i0+1)*n2+i1*n+i2];
                        T z2=Z[(i0+1)*n2+i1*n+i2];
                        line_segment_length[i0]=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
                        line_total_length+=line_segment_length[i0];
                    }
                    for (uint32_t i0=0;i0<n-1;i0++)
                    {
                        line_segment_length[i0]*=((T)2/line_total_length);
                        block_segment_length[i0]+=line_segment_length[i0];
                        block_segment_length_2[i0]+=line_segment_length[i0]*line_segment_length[i0];
                    }
                }
        }
        #pragma omp critical
        {
            for (uint32_t i=0;i<n-1;i++)
            {
                segment_length[i]+=block_segment_length[i];
                segment_length_2[i]+=block_segment_length_2[i];
            }
        }
    }
    vector<T> segment_mse(n-1,0);
    for (uint32_t i=0;i<n-1;i++)
    {
        segment_length[i]/=(params.element_num_to_process*n2*3);
        segment_length_2[i]/=(params.element_num_to_process*n2*3);
        segment_mse[i]=segment_length_2[i]-segment_length[i]*segment_length[i];
    }
    for (uint32_t i=0;i<n-1;i++)
        printf("%f ",segment_length[i]);
    printf("\n");
    for (uint32_t i=0;i<n-1;i++)
        printf("%f ",segment_mse[i]);
    printf("\n");
    exit(0);
}