#include "main.hpp"

template<typename T>
void compress1(jzzz::param_class<T> &params);
template<typename T>
void compress2(jzzz::param_class<T> &params,string field);

int main(int argc,char **argv)
{
    jzzz::param_class<float> params;
    for (int i=1;i<argc;i++)
    {
        string arg=argv[i];
        if (arg=="-i")
        {
            if (i+1<argc)
                params.data_path=argv[++i];
            else
                throw runtime_error("缺少-i参数后的输入文件路径");
            continue;
        }
        if (arg=="-elements")
        {
            if (i+1<argc)
            {
                string elements_str=argv[++i];
                elements_str=elements_str.substr(1,elements_str.size()-2);
                params.element_begin_index=stoull(elements_str.substr(0,elements_str.find(",")));
                params.element_end_index=stoull(elements_str.substr(elements_str.find(",")+1));
                params.got_element_range=true;
            }
            else
                throw runtime_error("缺少-elements参数后的element编号列表");
            continue;
        }
        if (arg=="-fields")
        {
            if (i+1<argc)
            {
                string fields_str=argv[++i];
                fields_str=fields_str.substr(1,fields_str.size()-2);
                stringstream ss(fields_str);
                string field_str;
                while (getline(ss,field_str,','))
                    params.fields_to_process.push_back(field_str);
            }
            else
                throw runtime_error("缺少-fields参数后的字段列表");
            continue;
        }
        if (arg=="-M")
        {
            if (i+2<argc)
            {
                params.eb_mode=argv[++i];
                params.input_eb=stof(argv[++i]);
            }
            else
                throw runtime_error("缺少-M参数后的误差控制模式和误差值");
            continue;
        }
        if (arg=="-method")
        {
            if (i+1<argc)
                params.method=argv[++i];
            else
                throw runtime_error("缺少-method参数后的压缩方法");
            continue;
        }
    }
    
    if (params.data_path.empty())
        throw runtime_error("必须指定原始数据路径，使用-i参数");
    if (params.eb_mode.empty())
        throw runtime_error("必须指定误差控制模式，使用-M参数");
    if (params.input_eb<=0)
        throw runtime_error("误差值必须为正数");
    if (params.fields_to_process.empty())
        params.all_fields_to_process=true;
    else
        params.all_fields_to_process=false;
    if (params.method.empty())
        params.method="jzzz";
    
    params.read_header(params.data_path);
    params.save_header_and_element_map_compress();
    params.get_max_cpus();

    params.blocks.resize(params.max_cpus);
    for (int i=0;i<params.max_cpus;i++)
    {
        auto &block=params.blocks[i];
        block.block_id=i;
        block.begin_index=(i*params.element_num_to_process)/params.max_cpus;
        block.end_index=((i+1)*params.element_num_to_process)/params.max_cpus;
        block.element_num=block.end_index-block.begin_index;
    }

    
    
    if (params.header.word_size==8)
    {
        jzzz::param_class<double> params_double(params);
        compress1(params_double);
    }
    else
        compress1(params);
}

template<typename T>
void compress1(jzzz::param_class<T> &params) //已经根据原始数据的格式改变了T
{
    for (auto field:params.fields_to_process)
    {
        printf("Field= %s\n",field.c_str());
        if (field=="X") //直接复制
        {
            auto header=params.header;
            ifstream data_file(params.data_path,ios::binary);
            data_file.seekg(params.field_begin_address[field],ios::beg);
            vector<T> coordinate_data(header.file_element_num*3*header.element_size);
            data_file.read(reinterpret_cast<char*>(coordinate_data.data()),header.file_element_num*3*header.element_size*header.word_size);
            data_file.close();
            ofstream coordinate_file(params.data_path+"_"+field+".jzzz",ios::binary);
            coordinate_file.write(reinterpret_cast<char*>(coordinate_data.data()),header.file_element_num*3*header.element_size*header.word_size);
            coordinate_file.close();
            continue;
        }
        params.data.resize(params.element_num_to_process*params.header.element_size);
        read_data(params,field);
        params.data_range_min[field]=params.data[0];
        params.data_range_max[field]=params.data[0];
        T global_min=params.data[0],global_max=params.data[0];
        #pragma omp parallel for schedule(static) num_threads(params.max_cpus) reduction(min:global_min) reduction(max:global_max)
        for (uint64_t block_id=0;block_id<params.blocks.size();block_id++)
        {
            auto &block=params.blocks[block_id];
            uint64_t begin_address=block.begin_index*params.header.element_size;
            uint64_t end_address=block.end_index*params.header.element_size;
            for (uint64_t i=begin_address;i<end_address;i++)
            {
                global_min=min(global_min,params.data[i]);
                global_max=max(global_max,params.data[i]);
            }
        }
        params.data_range_min[field]=global_min;
        params.data_range_max[field]=global_max;
        if (params.eb_mode=="ABS")
        {
            params.abs_eb[field]=params.input_eb;
            params.rel_eb[field]=params.input_eb/(params.data_range_max[field]-params.data_range_min[field]);
        }
        else if (params.eb_mode=="REL")
        {
            params.rel_eb[field]=params.input_eb;
            params.abs_eb[field]=params.input_eb*(params.data_range_max[field]-params.data_range_min[field]);
        }
        printf("Data range: [%f, %f], abs_eb=%e, rel_eb=%e\n",params.data_range_min[field],params.data_range_max[field],params.abs_eb[field],params.rel_eb[field]);

        if ((is_same<T,float>::value) and (params.rel_eb[field]<1e-6)) //原始数据类型是float，但是需要更高精度的压缩
        {
            jzzz::param_class<double> params_double(params);
            params_double.data.resize(params.element_num_to_process*params.header.element_size);
            #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
            for (uint64_t block_id=0;block_id<params.blocks.size();block_id++)
            {
                auto &block=params.blocks[block_id];
                uint64_t begin_address=block.begin_index*params.header.element_size;
                uint64_t end_address=block.end_index*params.header.element_size;
                for (uint64_t i=begin_address;i<end_address;i++)
                    params_double.data[i]=params.data[i];
            }
            params.data.clear();
            params.data.shrink_to_fit();
            compress2(params_double,field);
            params_double.data.clear();
            params_double.data.shrink_to_fit();
        }
        else
        {
            compress2(params,field);
            params.data.clear();
            params.data.shrink_to_fit();
        }
    }
}

template<typename T>
void compress2(jzzz::param_class<T> &params,string field) // 在最终决定params的T后，根据method调用不同的压缩函数，并写出压缩数据
{
    // get_coordinate_in_element(params); //用于统计每个element中相邻点的距离，验证了点确实在GLL点上。自带exit(0)，压缩时不使用。
    printf("params.method= %s\n",params.method.c_str());
    if (params.method=="jzzz")
        params.compress(field);
    else if (params.method=="sz2-0")
    {
        vector<T> backup_x=params.x;
        for (uint32_t i=0;i<params.x.size();i++)
            params.x[i]=(T)i/(params.x.size()-1)*2-1;
        sz2_1_compress(params,field);
        params.x=backup_x;
    }
    else if (params.method=="sz2-1")
        sz2_1_compress(params,field);
    else if (params.method=="sz3-1")
    {
        vector<T> backup_x=params.x;
        for (uint32_t i=0;i<params.x.size();i++)
            params.x[i]=(T)i/(params.x.size()-1)*2-1;
        sz3_2_compress(params,field);
        params.x=backup_x;
    }
    else if (params.method=="sz3-2")
        sz3_2_compress(params,field);
    else if (params.method=="zfp")
        zfp_compress(params,field);
    else
        throw runtime_error("未知的压缩方法: "+params.method);
        
    double br=0,cr=0;
    string compressed_path=params.data_path+"_"+field+".jzzz";
    ofstream compressed_file(compressed_path,ios::binary);
    uint8_t word_size_during_compression=sizeof(T);
    compressed_file.write(reinterpret_cast<char*>(&word_size_during_compression),1);
    uint64_t block_num=params.blocks.size();
    compressed_file.write(reinterpret_cast<char*>(&block_num),8);
    if (is_same<T,float>::value)
    {
        compressed_file.write(reinterpret_cast<char*>(&params.data_range_min[field]),4);
        compressed_file.write(reinterpret_cast<char*>(&params.data_range_max[field]),4);
    }
    else
    {
        float data_range_min_float=(float)(params.data_range_min[field]);
        float data_range_max_float=(float)(params.data_range_max[field]);
        compressed_file.write(reinterpret_cast<char*>(&data_range_min_float),4);
        compressed_file.write(reinterpret_cast<char*>(&data_range_max_float),4);
    }
    br+=17;
    
    for (auto &block:params.blocks)
    {
        uint64_t zstd_len=block.zstd_bs.size();
        compressed_file.write(reinterpret_cast<char*>(&zstd_len),8);
    }
    br+=8*params.blocks.size();

    for (auto &block:params.blocks)
    {
        compressed_file.write(reinterpret_cast<char*>(block.zstd_bs.data()),block.zstd_bs.size());
        br+=block.zstd_bs.size();
        block.zstd_bs.clear();
        block.zstd_bs.shrink_to_fit();
    }
    compressed_file.close();
    br*=((double)8/(params.element_num_to_process*params.header.element_size));
    cr=(8*params.header.word_size)/br;
    printf("\033[91mBR= %f bits/value, CR= %f\n\033[0m",br,cr);
}