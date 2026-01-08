#ifndef __READ_BLOCK_INFO_CPP__
#define __READ_BLOCK_INFO_CPP__

#include "jzzz_param.hpp"
#include <fstream>

namespace jzzz
{
    void param_class::read_block_info()
    {
        ifstream header_file(header_path);
        if (!header_file.is_open())
            throw runtime_error("无法打开header文件：" + header_path);

        string line;
        int line_index;
        for (line_index=0;;line_index++)
        {
            getline(header_file,line);
            if (line.find("blocks:")!=string::npos)
            {
                line_index++;
                break;
            }
        }
        
        bool file_end=false;
        for (;;line_index++)
        {
            if (!getline(header_file,line))
            {
                file_end=true;
                break;
            }
            block_class block;
            line.erase(0,string("- block_id: ").size());
            block.block_id=stoi(line);
            line_index++;
            getline(header_file,line);
            line.erase(0,string("  element_num: ").size());
            block.element_num=stoi(line);
            line_index++;
            
            for (;;line_index++)
            {
                if (!getline(header_file,line) || line.empty())
                {
                    file_end=true;
                    break;
                }
                if (line.find("- block_id: ")!=string::npos)
                {
                    header_file.seekg(-line.size()-1,ios::cur);
                    line_index--;
                    break;
                }
                line.erase(0,string("  ").size());
                string field=line.substr(0,line.find(": "));
                string path=line.substr(line.find(": ")+2);
                block.field_paths[field]=path;
            }
            this->blocks.push_back(block);
            if (file_end)
                break;
        }
        
        header_file.close();

        if (this->all_blocks_to_process)
            for (int i=0;i<int(this->blocks.size());i++)
                this->blocks_to_process.push_back(this->blocks[i].block_id);
        
        
    }
}

#endif