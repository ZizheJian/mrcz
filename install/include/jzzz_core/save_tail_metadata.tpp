#pragma once

namespace jzzz
{
    template<typename T>
    void param_class<T>::save_tail_metadata_compress()
    {
        ifstream data_file(data_path,ios::binary);
        uint64_t tail_begin=field_end_address[all_fields_in_file.back()];
        uint64_t tail_end;
        data_file.seekg(0,ios::end);
        tail_end=data_file.tellg();
        vector<unsigned char> tail_data(tail_end-tail_begin);
        data_file.seekg(tail_begin,ios::beg);
        data_file.read(reinterpret_cast<char*>(tail_data.data()),tail_end-tail_begin);
        data_file.close();

        ofstream tail_file(data_path+"_tail.jzzz",ios::binary);
        tail_file.write(reinterpret_cast<char*>(tail_data.data()),tail_data.size());
        tail_file.close();
    }

    template<typename T>
    void param_class<T>::save_tail_metadata_decompress()
    {
        ifstream tail_file(data_path+"_tail.jzzz",ios::binary);
        tail_file.seekg(0,ios::end);
        uint64_t tail_size=tail_file.tellg();
        tail_file.seekg(0,ios::beg);
        vector<unsigned char> tail_data(tail_size);
        tail_file.read(reinterpret_cast<char*>(tail_data.data()),tail_size);
        tail_file.close();

        ofstream decompressed_file(data_path+".jzzz.bin",ios::app|ios::binary);
        decompressed_file.seekp(0,ios::end);
        decompressed_file.write(reinterpret_cast<char*>(tail_data.data()),tail_size);
        decompressed_file.close();
    }
}