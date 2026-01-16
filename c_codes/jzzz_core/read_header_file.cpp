#ifndef __READ_HEADER_FILE_CPP__
#define __READ_HEADER_FILE_CPP__

#include "jzzz_param.hpp"
#include <fstream>

namespace jzzz
{
    void param_class::read_header_file()
    {
        ifstream header_file(this->header_path);
        if (!header_file.is_open())
            throw runtime_error("无法打开header文件：" + this->header_path);
        string line;
        
        for (int line_index=0;;line_index++)
        {
            getline(header_file,line);
            if (line.find("  word_size: ")!=string::npos)
            {
                line.erase(0,string("  word_size: ").size());
                this->header.word_size=stoi(line);
                continue;
            }
            if (line.find("  realtype: ")!=string::npos)
            {
                line.erase(0,string("  realtype: ").size());
                this->header.realtype=line;
                continue;
            }
            if (line.find("  element_shape: ")!=string::npos)
            {
                for (;;line_index++)
                {
                    getline(header_file,line);
                    if (line.find("  - ")!=string::npos)
                    {
                        line.erase(0,string("  - ").size());
                        this->header.element_shape.push_back(stoi(line));
                    }
                    else
                    {
                        header_file.seekg(-line.size()-1,ios::cur);
                        line_index--;
                        break;
                    }
                }
                continue;
            }
            if (line.find("  element_size: ")!=string::npos)
            {
                line.erase(0,string("  element_size: ").size());
                this->header.element_size=stoi(line);
                continue;
            }
            if (line.find("  file_element_num: ")!=string::npos)
            {
                line.erase(0,string("  file_element_num: ").size());
                this->header.file_element_num=stoi(line);;
                continue;
            }
            if (line.find("  variables: ")!=string::npos)
            {
                line.erase(0,string("  variables: ").size());
                this->header.variables=line;
                continue;
            }
            if (line.find("  passive_scalar_num: ")!=string::npos)
            {
                line.erase(0,string("  passive_scalar_num: ").size());
                this->header.passive_scalar_num=stoi(line);
                continue;
            }
            if (line.find("data_range:")!=string::npos)
                break;
        }

        header_file.close();

        if (this->all_fields_to_process)
            for (uint32_t i=0;i<this->header.variables.size();i++)
            {
                char c=this->header.variables[i];
                if (c=='U')
                {
                    this->fields_to_process.push_back("velocityx");
                    this->fields_to_process.push_back("velocityy");
                    this->fields_to_process.push_back("velocityz");
                }
                else if (c=='P')
                {
                    this->fields_to_process.push_back("pressure");
                }
                else if (c=='T')
                {
                    this->fields_to_process.push_back("temperature");
                }
                else if (c=='S')
                {
                    for (uint32_t j=0;j<this->header.passive_scalar_num;j++)
                        if (j<10)
                            this->fields_to_process.push_back("s0"+to_string(j));
                        else
                            this->fields_to_process.push_back("s"+to_string(j));
                    break;
                }
            }
    }
}

#endif