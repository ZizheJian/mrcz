#pragma once

#include "zfp_compress.tpp"

template<typename T>
void zfp_decompress_field(mrcz::param_class<T> &params,mrcz::partition_class &partition,string method,T *output,T abs_eb,const unsigned char *compressed,size_t compressed_size)
{
    zfp_stream *zfp=zfp_stream_open(nullptr);
    zfp_stream_set_accuracy(zfp,abs_eb);
    bitstream *input_bs=stream_open(const_cast<unsigned char*>(compressed),compressed_size);
    zfp_stream_set_bit_stream(zfp,input_bs);
    zfp_field *zfp_data=make_zfp_api_field(params,partition,method,output);
    if (!zfp_decompress(zfp,zfp_data))
        throw runtime_error("ZFP解压失败");
    stream_close(input_bs);
    zfp_field_free(zfp_data);
    zfp_stream_close(zfp);
}

template<typename T>
void zfp_decompress(mrcz::param_class<T> &params,string field)
{
    auto &header=params.header;
    auto &partitions=params.partitions;
    auto &decompressed_data=params.decompressed_data;
    T abs_eb=params.abs_eb[field];
    decompressed_data.resize(params.element_num_to_process*header.element_size);
    double start_time=omp_get_wtime();
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
    for (uint64_t partition_id=0;partition_id<partitions.size();partition_id++)
    {
        auto &partition=partitions[partition_id];
        uint64_t begin=partition.begin_index*header.element_size;
        if (params.method=="ZFP-3D-per-element")
        {
            size_t offset=0;
            for (uint64_t element_id=0;element_id<partition.element_num;element_id++)
            {
                uint64_t compressed_size=mrcz::baseline::read_uint64(partition.zstd_bs,offset);
                if (compressed_size>partition.zstd_bs.size()-offset)
                    throw runtime_error(params.method+"逐element bitstream不完整");
                zfp_decompress_field(params,partition,params.method,decompressed_data.data()+begin+element_id*header.element_size,abs_eb,partition.zstd_bs.data()+offset,compressed_size);
                offset+=compressed_size;
            }
            if (offset!=partition.zstd_bs.size())
                throw runtime_error(params.method+"逐element bitstream包含多余数据");
        }
        else
            zfp_decompress_field(params,partition,params.method,decompressed_data.data()+begin,abs_eb,partition.zstd_bs.data(),partition.zstd_bs.size());
    }
    double end_time=omp_get_wtime();
    printf("Decompression time= %f seconds\n",end_time-start_time);
    printf("\033[91mField= %s, decompression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),params.element_num_to_process*header.element_size*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
}
