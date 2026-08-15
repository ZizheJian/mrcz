#pragma once

namespace mrcz
{
    template<typename T>
    void param_class<T>::save_header_and_element_map(const string &task)
    {
        if (format=="nekrs")
            save_header_and_element_map_nekrs(task);
        else if (format=="mfem")
            save_header_and_element_map_mfem(task);
        else if (format=="nektar++")
            save_header_and_element_map_nektarpp(task);
    }

    template<typename T>
    void param_class<T>::save_header_and_element_map_nekrs(const string &task)
    {
        if (task=="compress")
        {
            ofstream header_file(data_path+"_header.mrcz");
            header_file.write(&header.original_header[0],header.original_header.size());
            header_file.close();
            string element_map;
            element_map.resize(header.file_element_num*4);
            ifstream data_file(data_path,ios::binary);
            data_file.seekg(136,ios::beg);
            data_file.read(&element_map[0],header.file_element_num*4);
            data_file.close();
            ofstream element_map_file(data_path+"_element_map.mrcz",ios::binary);
            element_map_file.write(&element_map[0],header.file_element_num*4);
            element_map_file.close();
        }
        else
        {
            ofstream decompressed_file(data_path+"."+method+".bin",ios::binary);
            decompressed_file.write(&header.original_header[0],header.original_header.size());
            string element_map;
            element_map.resize(header.file_element_num*4);
            ifstream element_map_file(data_path+"_element_map.mrcz",ios::binary);
            element_map_file.read(&element_map[0],header.file_element_num*4);
            element_map_file.close();
            decompressed_file.write(&element_map[0],header.file_element_num*4);
            decompressed_file.close();
        }
    }

    template<typename T>
    void param_class<T>::save_header_and_element_map_mfem(const string &task)
    {
        if (task=="compress")
        {
            // 只需要写header，没有element_map
            ofstream header_file(data_path+"_header.mrcz");
            header_file.write(header.original_header.data(),header.original_header.size());
            string s="DofNum="+to_string(header.dof_num)+"\n";
            header_file.write(s.data(),s.size());
            header_file.close();
        }
        else
        {
            // 什么都不干
        }
    }

    template<typename T>
    void param_class<T>::save_header_and_element_map_nektarpp(const string &task)
    {
        if (task=="compress")
            header.header_doc->SaveFile((data_path+"_header.mrcz").c_str());
        else
        {
            // 什么都不干
        }
    }
}
