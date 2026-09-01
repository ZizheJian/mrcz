#pragma once

#include "sperr_compress.tpp"

template<typename T>
void sperr_decompress_field(mrcz::param_class<T> &params,mrcz::partition_class &partition,string method,T *output,const unsigned char *compressed,size_t compressed_size)
{
    auto &header=params.header;
    size_t n0=header.c_element_shape[0];
    size_t n1=header.c_element_shape[1];
    size_t n2=header.c_element_shape[2];
    int is_float=is_same<T,float>::value?1:0;
    void *decompressed=nullptr;
    int status;
    if (method=="SPERR-2D")
        status=C_API::sperr_decomp_2d(compressed,compressed_size,is_float,header.element_size,partition.element_num,&decompressed);
    else if (method=="SPERR-3D-fused" or method=="SPERR-3D-per-element")
    {
        size_t dimx=0,dimy=0,dimz=0;
        status=C_API::sperr_decomp_3d(compressed,compressed_size,is_float,1,&dimx,&dimy,&dimz,&decompressed);
        size_t expected_element_num=method=="SPERR-3D-per-element"?1:partition.element_num;
        if (status==0 and (dimx!=n2 or dimy!=n1 or dimz!=expected_element_num*n0))
            throw runtime_error(params.method+" bitstream尺寸不匹配");
    }
    else
        throw runtime_error("未知的SPERR布局: "+method);
    if (status!=0)
        throw runtime_error("SPERR解压失败，状态码="+to_string(status));
    size_t value_num=method=="SPERR-3D-per-element"?header.element_size:partition.element_num*header.element_size;
    memcpy(output,decompressed,value_num*sizeof(T));
    free(decompressed);
}

template<typename T>
void sperr_decompress(mrcz::param_class<T> &params,string field)
{
    auto &header=params.header;
    auto &partitions=params.partitions;
    auto &decompressed_data=params.decompressed_data;
    decompressed_data.resize(params.element_num_to_process*header.element_size);
    double start_time=omp_get_wtime();
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
    for (uint64_t partition_id=0;partition_id<partitions.size();partition_id++)
    {
        auto &partition=partitions[partition_id];
        uint64_t begin=partition.begin_index*header.element_size;
        if (params.method=="SPERR-3D-per-element")
        {
            size_t offset=0;
            for (uint64_t element_id=0;element_id<partition.element_num;element_id++)
            {
                uint64_t element_size=mrcz::baseline::read_uint64(partition.zstd_bs,offset);
                if (element_size>partition.zstd_bs.size()-offset)
                    throw runtime_error(params.method+"逐element bitstream不完整");
                sperr_decompress_field(params,partition,params.method,decompressed_data.data()+begin+element_id*header.element_size,partition.zstd_bs.data()+offset,element_size);
                offset+=element_size;
            }
            if (offset!=partition.zstd_bs.size())
                throw runtime_error(params.method+"逐element bitstream包含多余数据");
        }
        else
            sperr_decompress_field(params,partition,params.method,decompressed_data.data()+begin,partition.zstd_bs.data(),partition.zstd_bs.size());
    }
    double end_time=omp_get_wtime();
    printf("Decompression time= %f seconds\n",end_time-start_time);
    printf("\033[91mField= %s, decompression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),params.element_num_to_process*header.element_size*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
}
