#ifndef __READ_DATA_RANGE_CPP__
#define __READ_DATA_RANGE_CPP__

#include "jzzz_param.hpp"
#include <fstream>

namespace jzzz
{
    void param_class::read_data_range()
    {
        ifstream header_file(header_path);
        if (!header_file.is_open())
            throw runtime_error("无法打开header文件：" + header_path);

        string line;
        int line_index;
        for (line_index=0;;line_index++)
        {
            getline(header_file,line);
            if (line.find("data_range:")!=string::npos)
            {
                line_index++;
                break;
            }
        }

        for (;;line_index++)
        {
            getline(header_file,line);
            if (line.find("block_num: ")!=string::npos)
                break;
            line.erase(0,string("  ").size());
            string variable=line.substr(0,line.find(":"));
            this->data_range[variable]=vector<vector<float>>();
            int dim_num=1;
            if (variable=="X" || variable=="U")
                dim_num=3;
            if (variable=="S")
                dim_num=this->header.passive_scalar_num;
            if (variable=="X" || variable=="U" || variable=="S")
            {
                for (int dim=0;dim<dim_num;dim++)
                {
                    getline(header_file,line);
                    line.erase(0,string("  - - ").size());
                    float min_val=stof(line);
                    line_index++;
                    getline(header_file,line);
                    line.erase(0,string("    - ").size());
                    float max_val=stof(line);
                    line_index++;
                    this->data_range[variable].push_back(vector<float>{min_val,max_val});
                }
                continue;
            }
            else
            {
                getline(header_file,line);
                line.erase(0,string("  - ").size());
                float min_val=stof(line);
                line_index++;
                getline(header_file,line);
                line.erase(0,string("  - ").size());
                float max_val=stof(line);
                line_index++;
                this->data_range[variable].push_back(vector<float>{min_val,max_val});
            }
            
        }

        header_file.close();
    }
}

#endif