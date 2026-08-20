#include "qoz_bridge.hpp"

#include "QoZ/api/sz.hpp"

using namespace std;

template<typename T>
void compress_4d(const T *data,size_t n0,size_t n1,size_t n2,size_t element_num,
    double abs_eb,vector<unsigned char> &compressed)
{
    // QoZ的最后一个维度变化最快，因此[element,n2,n1,n0]与MRCZ的内存顺序一致。
    QoZ::Config config(element_num,n2,n1,n0);
    config.errorBoundMode=QoZ::EB_ABS;
    config.absErrorBound=abs_eb;
    config.QoZ=4;
    config.openmp=false;

    size_t compressed_size=0;
    char *compressed_data=SZ_compress<T>(config,data,compressed_size);
    compressed.assign(compressed_data,compressed_data+compressed_size);
    delete[] compressed_data;
}

template<typename T>
void decompress_4d(const vector<unsigned char> &compressed,T *data)
{
    QoZ::Config config;
    T *output=data;
    SZ_decompress<T>(config,reinterpret_cast<char *>(const_cast<unsigned char *>(compressed.data())),
        compressed.size(),output);
}

void qoz_compress_4d(const float *data,size_t n0,size_t n1,size_t n2,size_t element_num,
    double abs_eb,vector<unsigned char> &compressed)
{
    compress_4d(data,n0,n1,n2,element_num,abs_eb,compressed);
}

void qoz_compress_4d(const double *data,size_t n0,size_t n1,size_t n2,size_t element_num,
    double abs_eb,vector<unsigned char> &compressed)
{
    compress_4d(data,n0,n1,n2,element_num,abs_eb,compressed);
}

void qoz_decompress_4d(const vector<unsigned char> &compressed,float *data)
{
    decompress_4d(compressed,data);
}

void qoz_decompress_4d(const vector<unsigned char> &compressed,double *data)
{
    decompress_4d(compressed,data);
}
