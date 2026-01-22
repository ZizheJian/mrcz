#pragma once

namespace jzzz
{
    template<typename T>
    void param_class<T>::read_header(string header_path)
    {
        //////// 读取header文件 ////////
        ifstream data_file(header_path,ios::binary);
        header.original_header.resize(136);
        data_file.read(&header.original_header[0],136);
        data_file.close();
        //////// 检查是大端数据还是小端数据 ////////
        const char *etagb=header.original_header.data()+132;
        float etagL=0;
        memcpy(&etagL,etagb,4);
        etagL=static_cast<float>(static_cast<int>(etagL*1e5))/1e5;
        char etagb_reversed[4]={etagb[3],etagb[2],etagb[1],etagb[0]};
        float etagB=0;
        memcpy(&etagB,etagb_reversed,4);
        etagB=static_cast<float>(static_cast<int>(etagB*1e5))/1e5;
        if (fabs(etagL-6.54321)<1e-5)
            header.endian_mode='<';
        else if (fabs(etagB-6.54321)<1e-5)
            header.endian_mode='>';
        else
            throw runtime_error("无法识别数据是大端还是小端格式");
        //////// 读取剩余header信息 ////////
        istringstream header_stream(header.original_header);
        string word;
        int word_index=0;
        while (header_stream>>word)
        {
            if (word_index==1)
                header.word_size=stoi(word);
            else if (word_index==2)
                header.element_shape.push_back(stoi(word));
            else if (word_index==3)
                header.element_shape.push_back(stoi(word));
            else if (word_index==4)
                header.element_shape.push_back(stoi(word));
            else if (word_index==6)
                header.file_element_num=stoull(word);
            else if (word_index==11)
                header.variables=word;
            word_index++;
        }
        header.element_size=header.element_shape[0]*header.element_shape[1]*header.element_shape[2];
        size_t s_pos=header.variables.find('S');
        if (s_pos!=string::npos)
            header.passive_scalar_num=stoi(header.variables.substr(s_pos+1,2));
        else
            header.passive_scalar_num=0;

        //////// 生成all_fields_in_file,fields_to_process，并计算每个field在源文件中的起始和结束地址 ////////
        uint64_t last_field_end_address=136+header.file_element_num*4;
        for (uint32_t i=0;i<header.variables.size();i++)
        {
            char c=header.variables[i];
            if (c=='X')
            {
                all_fields_in_file.push_back("X");
                field_begin_address["X"]=last_field_end_address;
                last_field_end_address+=header.file_element_num*header.element_size*header.word_size*3;
                field_end_address["X"]=last_field_end_address;
            }
            else if (c=='U')
            {
                all_fields_in_file.push_back("UX");
                all_fields_in_file.push_back("UY");
                all_fields_in_file.push_back("UZ");
                field_begin_address["UX"]=last_field_end_address;
                field_begin_address["UY"]=last_field_end_address;
                field_begin_address["UZ"]=last_field_end_address;
                last_field_end_address+=header.file_element_num*header.element_size*header.word_size*3;
                field_end_address["UX"]=last_field_end_address;
                field_end_address["UY"]=last_field_end_address;
                field_end_address["UZ"]=last_field_end_address;
            }
            else if (c=='P')
            {
                all_fields_in_file.push_back("P");
                field_begin_address["P"]=last_field_end_address;
                last_field_end_address+=header.file_element_num*header.element_size*header.word_size;
                field_end_address["P"]=last_field_end_address;
            }
            else if (c=='T')
            {
                all_fields_in_file.push_back("T");
                field_begin_address["T"]=last_field_end_address;
                last_field_end_address+=header.file_element_num*header.element_size*header.word_size;
                field_end_address["T"]=last_field_end_address;
            }
            else if (c=='S')
            {
                for (uint32_t j=0;j<header.passive_scalar_num;j++)
                {
                    if (j<10)
                        all_fields_in_file.push_back("S0"+to_string(j));
                    else
                        all_fields_in_file.push_back("S"+to_string(j));
                    field_begin_address[all_fields_in_file.back()]=last_field_end_address;
                    last_field_end_address+=header.file_element_num*header.element_size*header.word_size;
                    field_end_address[all_fields_in_file.back()]=last_field_end_address;
                }
                i+=2;
            }
        }
        if (all_fields_to_process)
            fields_to_process=all_fields_in_file;
        
        //////// 计算实际需要处理的element范围 ////////
        if (!got_element_range)
        {
            element_begin_index=0;
            element_end_index=header.file_element_num;
            element_num_to_process=header.file_element_num;
        }
        else
            element_num_to_process=element_end_index-element_begin_index;
        
        //////// 调整fields的起始和结束地址 ////////
        for (auto field:fields_to_process)
            if (field=="UX" || field=="UY" || field=="UZ")
            {
                field_begin_address[field]+=header.element_size*header.word_size*3*element_begin_index;
                field_end_address[field]=field_begin_address[field]+header.element_size*header.word_size*3*element_num_to_process;
            }
            else
            {
                field_begin_address[field]+=header.element_size*header.word_size*element_begin_index;
                field_end_address[field]=field_begin_address[field]+header.element_size*header.word_size*element_num_to_process;
            }
    }
}