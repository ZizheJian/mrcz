#ifndef __JZZZ_PARAM_H__
#define __JZZZ_PARAM_H__

#include <string>
#include <vector>
#include <unordered_map>
using namespace std;

namespace jzzz
{
    class header_class
    {
        public:
            int word_size;
            string realtype;
            vector<int> element_shape;
            int element_size;
            int file_element_num;
            string variables;
            int passive_scalar_num;
    };

    class block_class
    {
        public:
            int block_id;
            int element_num;
            unordered_map<string,string> field_paths;
    };

    class param_class
    {
        public:
            string header_path;
            string eb_mode;
            float eb;
            bool all_blocks_to_process=true;
            vector<int> blocks_to_process;
            bool all_fields_to_process=true;
            vector<string> fields_to_process;

            header_class header;
            unordered_map<string,vector<vector<float>>> data_range;
            int block_num;
            vector<block_class> blocks;

            vector<float> data;

            void read_header_file();
            void read_data_range();
            void read_block_info();
            string get_full_field_name(string field);
            
    };
}

#endif