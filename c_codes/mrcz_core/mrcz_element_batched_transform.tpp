#pragma once

namespace mrcz
{
    template<typename T>
    inline void param_class<T>::forward_transform_element_batched(vector<T> &temp_element1,vector<T> &temp_element2,T *element_batch,uint32_t batch_width,uint32_t n0,uint32_t n1,uint32_t n2,uint32_t n12,uint32_t n012,
        vector<T> &L0,vector<T> &L1,vector<T> &L2)
    {
        T *temp1=temp_element1.data(),*temp2=temp_element2.data();
        T *l0=L0.data(),*l1=L1.data(),*l2=L2.data();
        memset(temp1,0,sizeof(T)*n012*batch_width);
        for (uint32_t a=0;a<n0;a++)
            for (uint32_t i=0;i<n0;i++)
            {
                T weight=l0[a*n0+i];
                for (uint32_t bc=0;bc<n12;bc++)
                {
                    T *input=element_batch+(a*n12+bc)*batch_width;
                    T *output=temp1+(i*n12+bc)*batch_width;
                    #pragma omp simd
                    for (uint32_t lane=0;lane<batch_width;lane++)
                        output[lane]+=input[lane]*weight;
                }
            }
        memset(temp2,0,sizeof(T)*n012*batch_width);
        for (uint32_t a=0;a<n0;a++)
            for (uint32_t b=0;b<n1;b++)
                for (uint32_t j=0;j<n1;j++)
                {
                    T weight=l1[b*n1+j];
                    for (uint32_t c=0;c<n2;c++)
                    {
                        T *input=temp1+((a*n1+b)*n2+c)*batch_width;
                        T *output=temp2+((a*n1+j)*n2+c)*batch_width;
                        #pragma omp simd
                        for (uint32_t lane=0;lane<batch_width;lane++)
                            output[lane]+=input[lane]*weight;
                    }
                }
        memset(temp1,0,sizeof(T)*n012*batch_width);
        for (uint32_t ab=0;ab<n0*n1;ab++)
            for (uint32_t c=0;c<n2;c++)
                for (uint32_t k=0;k<n2;k++)
                {
                    T weight=l2[c*n2+k];
                    T *input=temp2+(ab*n2+c)*batch_width;
                    T *output=temp1+(ab*n2+k)*batch_width;
                    #pragma omp simd
                    for (uint32_t lane=0;lane<batch_width;lane++)
                        output[lane]+=input[lane]*weight;
                }
    }

    template<typename T>
    inline void param_class<T>::inverse_transform_element_batched(vector<T> &temp_element1,vector<T> &temp_element2,uint32_t batch_width,uint32_t n0,uint32_t n1,uint32_t n2,uint32_t n12,uint32_t n012,
        vector<T> &neg_L0,vector<T> &neg_L1,vector<T> &neg_L2)
    {
        T *temp1=temp_element1.data(),*temp2=temp_element2.data();
        T *l0=neg_L0.data(),*l1=neg_L1.data(),*l2=neg_L2.data();
        memset(temp2,0,sizeof(T)*n012*batch_width);
        for (uint32_t ab=0;ab<n0*n1;ab++)
            for (uint32_t k=0;k<n2;k++)
                for (uint32_t c=0;c<n2;c++)
                {
                    T weight=l2[k*n2+c];
                    T *input=temp1+(ab*n2+k)*batch_width;
                    T *output=temp2+(ab*n2+c)*batch_width;
                    #pragma omp simd
                    for (uint32_t lane=0;lane<batch_width;lane++)
                        output[lane]+=input[lane]*weight;
                }
        memset(temp1,0,sizeof(T)*n012*batch_width);
        for (uint32_t a=0;a<n0;a++)
            for (uint32_t j=0;j<n1;j++)
                for (uint32_t b=0;b<n1;b++)
                {
                    T weight=l1[j*n1+b];
                    for (uint32_t c=0;c<n2;c++)
                    {
                        T *input=temp2+((a*n1+j)*n2+c)*batch_width;
                        T *output=temp1+((a*n1+b)*n2+c)*batch_width;
                        #pragma omp simd
                        for (uint32_t lane=0;lane<batch_width;lane++)
                            output[lane]+=input[lane]*weight;
                    }
                }
        memset(temp2,0,sizeof(T)*n012*batch_width);
        for (uint32_t i=0;i<n0;i++)
            for (uint32_t a=0;a<n0;a++)
            {
                T weight=l0[i*n0+a];
                for (uint32_t bc=0;bc<n12;bc++)
                {
                    T *input=temp1+(i*n12+bc)*batch_width;
                    T *output=temp2+(a*n12+bc)*batch_width;
                    #pragma omp simd
                    for (uint32_t lane=0;lane<batch_width;lane++)
                        output[lane]+=input[lane]*weight;
                }
            }
    }
}
