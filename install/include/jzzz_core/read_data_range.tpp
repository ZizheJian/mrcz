#pragma once

namespace jzzz
{
    template<typename T>
    void param_class<T>::read_data_range()
    {
        // ifstream header_file(yaml_path);
        // if (!header_file.is_open())
        //     throw runtime_error("无法打开header文件：" + yaml_path);

        // string line;
        // int line_index;
        // for (line_index=0;;line_index++)
        // {
        //     getline(header_file,line);
        //     if (line.find("data_range:")!=string::npos)
        //     {
        //         line_index++;
        //         break;
        //     }
        // }

        // for (;;line_index++)
        // {
        //     getline(header_file,line);
        //     if (line.find("block_num: ")!=string::npos)
        //         break;
        //     line.erase(0,string("  ").size());
        //     string variable=line.substr(0,line.find(":"));
        //     this->data_range[variable]=vector<vector<T>>();
        //     int dim_num=1;
        //     if (variable=="X" || variable=="U")
        //         dim_num=3;
        //     if (variable=="S")
        //         dim_num=this->header.passive_scalar_num;
        //     if (variable=="X" || variable=="U" || variable=="S")
        //     {
        //         for (int dim=0;dim<dim_num;dim++)
        //         {
        //             getline(header_file,line);
        //             line.erase(0,string("  - - ").size());
        //             T min_val=stof(line);
        //             line_index++;
        //             getline(header_file,line);
        //             line.erase(0,string("    - ").size());
        //             T max_val=stof(line);
        //             line_index++;
        //             this->data_range[variable].push_back(vector<T>{min_val,max_val});
        //         }
        //         continue;
        //     }
        //     else
        //     {
        //         getline(header_file,line);
        //         line.erase(0,string("  - ").size());
        //         T min_val=stof(line);
        //         line_index++;
        //         getline(header_file,line);
        //         line.erase(0,string("  - ").size());
        //         T max_val=stof(line);
        //         line_index++;
        //         this->data_range[variable].push_back(vector<T>{min_val,max_val});
        //     }
            
        // }

        // header_file.close();
    }
}