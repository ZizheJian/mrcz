#pragma once

#include "mrcz_core/mrcz_param.hpp"

template<typename T>
void read_data(mrcz::param_class<T> &params,string field)//分配好params.data的空间后，调用这个函数来读取指定field的数据到params.data中
{
    auto &header=params.header;
    ifstream data_file(params.data_path,ios::binary);
    data_file.seekg(params.field_begin_address[field],ios::beg);
    if (field=="UX" or field=="UY" or field=="UZ")
    {
        vector<T> temp_data(params.element_num_to_process*header.element_size*3);
        data_file.read(reinterpret_cast<char*>(temp_data.data()),params.element_num_to_process*header.element_size*3*header.word_size);
        if (field=="UX")
            for (uint64_t element_id=0;element_id<params.element_num_to_process;element_id++)
                memcpy(&params.data[element_id*header.element_size],&temp_data[3*element_id*header.element_size],
                    header.element_size*header.word_size);
        if (field=="UY")
            for (uint64_t element_id=0;element_id<params.element_num_to_process;element_id++)
                memcpy(&params.data[element_id*header.element_size],&temp_data[(3*element_id+1)*header.element_size],
                    header.element_size*header.word_size);
        if (field=="UZ")
            for (uint64_t element_id=0;element_id<params.element_num_to_process;element_id++)
                memcpy(&params.data[element_id*header.element_size],&temp_data[(3*element_id+2)*header.element_size],
                    header.element_size*header.word_size);
    }
    else
        data_file.read(reinterpret_cast<char*>(params.data.data()),params.element_num_to_process*header.element_size*header.word_size);
    data_file.close();
}

template<typename T>
void write_data(mrcz::param_class<T> &params,string field)
{
    assert(params.header.word_size==sizeof(T));
    fstream decompressed_file(params.data_path+".mrcz.bin",ios::app|ios::binary);
    if (field=="UX" || field=="UY" || field=="UZ")
    {
        if (field=="UX")
        {
            decompressed_file.seekp(0,ios::end);
            assert(decompressed_file.tellp()==params.field_begin_address[field]);
            vector<T> temp_data(params.element_num_to_process*params.header.element_size*3,0);
            for (uint64_t element_id=0;element_id<params.element_num_to_process;element_id++)
                memcpy(&temp_data[element_id*params.header.element_size*3],&params.decompressed_data[element_id*params.header.element_size],params.header.element_size*sizeof(T));
            decompressed_file.write(reinterpret_cast<char*>(temp_data.data()),params.element_num_to_process*params.header.element_size*3*sizeof(T));
        }
        else if (field=="UY")
        {
            decompressed_file.seekp(0,ios::end);
            assert(decompressed_file.tellp()==params.field_end_address[field]);
            vector<T> temp_data(params.element_num_to_process*params.header.element_size*3,0);
            decompressed_file.seekg(params.field_begin_address[field]);
            decompressed_file.read(reinterpret_cast<char*>(temp_data.data()),params.element_num_to_process*params.header.element_size*3*sizeof(T));
            for (uint64_t element_id=0;element_id<params.element_num_to_process;element_id++)
                memcpy(&temp_data[(3*element_id+1)*params.header.element_size],&params.decompressed_data[element_id*params.header.element_size],
                    params.header.element_size*sizeof(T));
            decompressed_file.write(reinterpret_cast<char*>(temp_data.data()),params.element_num_to_process*params.header.element_size*3*sizeof(T));
        }
        else if (field=="UZ")
        {
            decompressed_file.seekp(0,ios::end);
            assert(decompressed_file.tellp()==params.field_end_address[field]);
            vector<T> temp_data(params.element_num_to_process*params.header.element_size*3,0);
            decompressed_file.seekg(params.field_begin_address[field]);
            decompressed_file.read(reinterpret_cast<char*>(temp_data.data()),params.element_num_to_process*params.header.element_size*3*sizeof(T));
            for (uint64_t element_id=0;element_id<params.element_num_to_process;element_id++)
                memcpy(&temp_data[(3*element_id+2)*params.header.element_size],&params.decompressed_data[element_id*params.header.element_size],
                    params.header.element_size*sizeof(T));
            decompressed_file.write(reinterpret_cast<char*>(temp_data.data()),params.element_num_to_process*params.header.element_size*3*sizeof(T));
        }
    }
    else
    {
        decompressed_file.seekp(0,ios::end);
        assert(decompressed_file.tellp()==params.field_begin_address[field]);
        decompressed_file.write(reinterpret_cast<char*>(params.decompressed_data.data()),params.element_num_to_process*params.header.element_size*params.header.word_size);
    }
    decompressed_file.close();
}