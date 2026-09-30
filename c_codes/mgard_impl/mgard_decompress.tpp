#pragma once

#include "mgard_compress.tpp"

template<typename T>
void mgard_decompress_array(const unsigned char *compressed,size_t compressed_size,T *output,const vector<mgard_x::SIZE> &expected_shape)
{
    void *output_ptr=output;
    vector<mgard_x::SIZE> shape;
    mgard_x::data_type dtype;
    auto status=mgard_x::decompress(compressed,compressed_size,output_ptr,shape,dtype,make_mgard_config(),true);
    if(status!=mgard_x::compress_status_type::Success)
        throw runtime_error("MGARD解压失败，状态码="+to_string(static_cast<int>(status)));
    auto expected_dtype=is_same_v<T,float>?mgard_x::data_type::Float:mgard_x::data_type::Double;
    if(output_ptr!=output or shape!=expected_shape or dtype!=expected_dtype)
        throw runtime_error("MGARD解压结果的地址、尺寸或数据类型不匹配");
}

template<typename T>
void mgard_decompress(mrcz::param_class<T> &params,string field)
{
    auto &header=params.header;
    auto &partitions=params.partitions;
    auto &decompressed_data=params.decompressed_data;
    decompressed_data.resize(params.element_num_to_process*header.element_size);
    double start_time=omp_get_wtime();
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
    for(uint64_t partition_id=0;partition_id<partitions.size();partition_id++)
    {
        auto &partition=partitions[partition_id];
        uint64_t begin=partition.begin_index*header.element_size;
        auto grid=make_mgard_api_grid(params,partition,params.method);
        if(partition.zstd_bs.empty())
            throw runtime_error(params.method+" partition bitstream为空");
        if(mgard_is_per_element(params.method))
        {
            size_t offset=0;
            for(uint64_t element_id=0;element_id<partition.element_num;element_id++)
            {
                uint64_t compressed_size=mrcz::baseline::read_uint64(partition.zstd_bs,offset);
                if(compressed_size>partition.zstd_bs.size()-offset)
                    throw runtime_error(params.method+"逐element bitstream不完整");
                mgard_decompress_array(partition.zstd_bs.data()+offset,compressed_size,decompressed_data.data()+begin+element_id*header.element_size,grid.shape);
                offset+=compressed_size;
            }
            if(offset!=partition.zstd_bs.size())
                throw runtime_error(params.method+"逐element bitstream包含多余数据");
        }
        else
            mgard_decompress_array(partition.zstd_bs.data(),partition.zstd_bs.size(),decompressed_data.data()+begin,grid.shape);
    }
    double end_time=omp_get_wtime();
    printf("Decompression time= %f seconds\n",end_time-start_time);
    printf("\033[91mField= %s, decompression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),params.element_num_to_process*header.element_size*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
}
