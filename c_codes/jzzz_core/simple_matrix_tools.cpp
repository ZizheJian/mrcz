#ifndef __SIMPLE_MATRIX_TOOLS_CPP__
#define __SIMPLE_MATRIX_TOOLS_CPP__

#include "jzzz_param.hpp"
#include <stdexcept>

namespace jzzz
{
    vector<float> param_class::simple_invert_matrix(vector<float>& A,uint32_t n)
    {
        vector<float> mat(2*n*n,0);
        for (uint32_t i=0;i<n;i++)
        {
            for (uint32_t j=0;j<n;j++)
            {
                mat[i*2*n+j]=A[i*n+j];
            }
            mat[i*2*n+i+n]=1;
        }
        for (uint32_t p=0;p<n;p++)
        {
            uint32_t max_row=p;
            for (uint32_t i=p+1;i<n;++i)
                if (abs(mat[i*2*n+p])>abs(mat[max_row*2*n+p]))
                    max_row=i;
            if (p!=max_row)
                for (uint32_t i=p;i<2*n;i++)
                    swap(mat[p*2*n+i],mat[max_row*2*n+i]);
            if (abs(mat[p*2*n+p])<1e-12)
                throw runtime_error("Matrix is singular or nearly singular");
            // Normalize pivot row
            float pivot=mat[p*2*n+p];
            for (uint32_t j=p;j<2*n;j++)
                mat[p*2*n+j]/=pivot;
            // Eliminate
            for (uint32_t i=0;i<n;i++)
            {
                if (i==p)
                    continue;
                float factor=mat[i*2*n+p];
                for (uint32_t j=p;j<2*n;++j)
                    mat[i*2*n+j]-=factor*mat[p*2*n+j];
            }
        }
        // Extract inverse
        vector<float> inv(n*n,0.0);
        for (uint32_t i=0;i<n;++i)
            for (uint32_t j=0;j<n;++j)
                inv[i*n+j]=mat[i*2*n+j+n];
        return inv;
    }

    vector<float> param_class::simple_matrix_multiply(vector<float>& A,vector<float>& B,uint32_t n)
    {
        vector<float> C(n*n,0.0);
        for (uint32_t i=0;i<n;i++)
            for (uint32_t j=0;j<n;j++)
                for (uint32_t k=0;k<n;k++)
                    C[i*n+j]+=A[i*n+k]*B[k*n+j];
        return C;
    }

}

#endif