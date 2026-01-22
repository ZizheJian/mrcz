#pragma once

namespace jzzz
{
    template<typename T>
    void param_class<T>::save_header_and_element_map_compress()
    {
        ofstream header_file(data_path+"_header.jzzz");
        header_file.write(&header.original_header[0],header.original_header.size());
        header_file.close();
        string element_map;
        element_map.resize(header.file_element_num*4);
        ifstream data_file(data_path,ios::binary);
        data_file.seekg(136,ios::beg);
        data_file.read(&element_map[0],header.file_element_num*4);
        data_file.close();
        ofstream element_map_file(data_path+"_element_map.jzzz",ios::binary);
        element_map_file.write(&element_map[0],header.file_element_num*4);
        element_map_file.close();
    }

    template<typename T>
    void param_class<T>::save_header_and_element_map_decompress()
    {
        ofstream decompressed_file(data_path+".jzzz.bin",ios::binary);
        decompressed_file.write(&header.original_header[0],header.original_header.size());
        string element_map;
        element_map.resize(header.file_element_num*4);
        ifstream element_map_file(data_path+"_element_map.jzzz",ios::binary);
        element_map_file.read(&element_map[0],header.file_element_num*4);
        element_map_file.close();
        decompressed_file.write(&element_map[0],header.file_element_num*4);
        decompressed_file.close();
    }
}