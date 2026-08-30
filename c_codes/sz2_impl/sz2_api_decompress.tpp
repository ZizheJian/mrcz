#pragma once

#include "sz2_api_compress.tpp"

template<typename T>
void sz2_decompress_field(mrcz::param_class<T> &params,mrcz::partition_class &partition,mrcz::baseline::api_layout layout,T *output,const unsigned char *compressed,size_t compressed_size)
{
    sz2_api_dims dims=make_sz2_api_dims(params,partition,layout);
    int data_type=is_same<T,float>::value?SZ_FLOAT:SZ_DOUBLE;
    void *decompressed=SZ_decompress(data_type,const_cast<unsigned char*>(compressed),compressed_size,dims.r5,dims.r4,dims.r3,dims.r2,dims.r1);
    if (decompressed==nullptr)
        throw runtime_error("SZ2解压失败");
    size_t value_num=layout==mrcz::baseline::api_layout::per_element_3d?params.header.element_size:partition.element_num*params.header.element_size;
    memcpy(output,decompressed,value_num*sizeof(T));
    free(decompressed);
}

template<typename T>
void sz2_api_decompress(mrcz::param_class<T> &params,string field)
{
    auto &header=params.header;
    auto &partitions=params.partitions;
    auto &decompressed_data=params.decompressed_data;
    auto layout=mrcz::baseline::parse_api_layout(params.method,"SZ2");
    decompressed_data.resize(params.element_num_to_process*header.element_size);
    if (SZ_Init(nullptr)!=SZ_SCES)
        throw runtime_error("SZ2初始化失败");
    double start_time=omp_get_wtime();
    // XXX SZ2官方threadsafe解压入口仍直接调用共享全局状态的普通API，因此必须串行处理partition
    for (uint64_t partition_id=0;partition_id<partitions.size();partition_id++)
    {
        auto &partition=partitions[partition_id];
        uint64_t begin=partition.begin_index*header.element_size;
        if (layout==mrcz::baseline::api_layout::per_element_3d)
        {
            size_t offset=0;
            uint64_t element_num=mrcz::baseline::read_uint64(partition.zstd_bs,offset);
            if (element_num!=partition.element_num)
                throw runtime_error(params.method+"逐element bitstream中的element数量不匹配");
            for (uint64_t element_id=0;element_id<element_num;element_id++)
            {
                uint64_t element_size=mrcz::baseline::read_uint64(partition.zstd_bs,offset);
                if (element_size>partition.zstd_bs.size()-offset)
                    throw runtime_error(params.method+"逐element bitstream不完整");
                sz2_decompress_field(params,partition,layout,decompressed_data.data()+begin+element_id*header.element_size,partition.zstd_bs.data()+offset,element_size);
                offset+=element_size;
            }
            if (offset!=partition.zstd_bs.size())
                throw runtime_error(params.method+"逐element bitstream包含多余数据");
        }
        else
            sz2_decompress_field(params,partition,layout,decompressed_data.data()+begin,partition.zstd_bs.data(),partition.zstd_bs.size());
    }
    double end_time=omp_get_wtime();
    SZ_Finalize();
    printf("Decompression time= %f seconds\n",end_time-start_time);
    printf("\033[91mField= %s, decompression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),params.element_num_to_process*header.element_size*header.word_size/((end_time-start_time)*1024*1024));
}
