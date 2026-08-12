#include "main.hpp"

template<typename T>
void compress1(mrcz::param_class<T> &params);
template<typename T>
void compress2(mrcz::param_class<T> &params,string field);

int main(int argc,char **argv)
{
    mrcz::param_class<float> params;
    for (int i=1;i<argc;i++)
    {
        string arg=argv[i];
        if (arg=="-format")
        {
            if (i+1<argc)
                params.format=argv[++i];
            else
                throw runtime_error("缺少-format参数后的数据格式");
            continue;
        }
        if (arg=="-i")
        {
            if (i+1<argc)
                params.data_path=argv[++i];
            else
                throw runtime_error("缺少-i参数后的输入文件路径");
            continue;
        }
        if (arg=="-mesh")
        {
            if (i+1<argc)
                params.mesh_path=argv[++i];
            else
                throw runtime_error("缺少-mesh参数后的网格文件路径");
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
                throw runtime_error("缺少-fields参数后的物理量列表");
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
    
    if (params.format.empty() or (params.format!="nekrs" and params.format!="mfem" and params.format!="nektar++"))
        throw runtime_error("必须指定数据格式，使用-format参数，目前支持的格式有: nekrs, mfem, nektar++");
    if (params.data_path.empty())
        throw runtime_error("必须指定原始数据路径，使用-i参数");
    if (params.format=="mfem" and params.mesh_path.empty())
        throw runtime_error("使用mfem格式时必须指定网格文件路径，使用-mesh参数");
    if (params.eb_mode.empty())
        throw runtime_error("必须指定误差控制模式，使用-M参数");
    if (params.input_eb<=0)
        throw runtime_error("误差值必须为正数");
    if (params.fields_to_process.empty())
        params.all_fields_to_process=true;
    else
        params.all_fields_to_process=false;
    if (params.method.empty())
        params.method="mrcz";

    params.read_header("compress");
    params.default_init();
    params.save_header_and_element_map("compress");
    params.max_cpus=get_max_cpus(params.max_cpus);

    params.partitions.resize(params.max_cpus);
    for (int i=0;i<params.max_cpus;i++)
    {
        auto &partition=params.partitions[i];
        partition.partition_id=i;
        partition.begin_index=(i*params.element_num_to_process)/params.max_cpus;
        partition.end_index=((i+1)*params.element_num_to_process)/params.max_cpus;
        partition.element_num=partition.end_index-partition.begin_index;
    }
    
    if (params.header.word_size==8)
    {
        mrcz::param_class<double> params_double(params);
        compress1(params_double);
    }
    else
        compress1(params);
}

template<typename T>
void compress1(mrcz::param_class<T> &params) //已经根据原始数据的格式改变了T
{
    for (auto field:params.fields_to_process)
    {
        printf("Field= %s\n",field.c_str());
        // if (field=="X") //直接复制
        // {
        //     auto header=params.header;
        //     ifstream data_file(params.data_path,ios::binary);
        //     data_file.seekg(params.field_begin_address[field],ios::beg);
        //     vector<T> coordinate_data(header.file_element_num*3*header.element_size);
        //     data_file.read(reinterpret_cast<char*>(coordinate_data.data()),header.file_element_num*3*header.element_size*header.word_size);
        //     data_file.close();
        //     ofstream coordinate_file(params.data_path+"_"+field+".mrcz",ios::binary);
        //     coordinate_file.write(reinterpret_cast<char*>(coordinate_data.data()),header.file_element_num*3*header.element_size*header.word_size);
        //     coordinate_file.close();
        //     continue;
        // }
        params.data.resize(params.element_num_to_process*params.header.element_size);
        read_data(params,field);
        T global_min=params.data[0],global_max=params.data[0];
        #pragma omp parallel for schedule(static) num_threads(params.max_cpus) reduction(min:global_min) reduction(max:global_max)
        for (uint64_t partition_id=0;partition_id<params.partitions.size();partition_id++)
        {
            auto &partition=params.partitions[partition_id];
            uint64_t begin_address=partition.begin_index*params.header.element_size;
            uint64_t end_address=partition.end_index*params.header.element_size;
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
        
        if ((is_same<T,float>::value) and (params.rel_eb[field]<1e-4)) //原始数据类型是float，但是需要更高精度的压缩
        {
            mrcz::param_class<double> params_double(params);
            params_double.data.resize(params.element_num_to_process*params.header.element_size);
            #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
            for (uint64_t partition_id=0;partition_id<params.partitions.size();partition_id++)
            {
                auto &partition=params.partitions[partition_id];
                uint64_t begin_address=partition.begin_index*params.header.element_size;
                uint64_t end_address=partition.end_index*params.header.element_size;
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
void compress2(mrcz::param_class<T> &params,string field) // 在最终决定params的T后，根据method调用不同的压缩函数，并写出压缩数据
{
    // tools/get_coordinate_in_element.tpp用于统计每个element中相邻点的距离，验证了点确实在GLL点上。自带exit(0)，压缩时不使用。
    if (params.data_range_min==params.data_range_max)
        printf("Warning: field %s has constant value %f\n",field.c_str(),params.data_range_min[field]);
    else
    {
        vector<vector<T>> backup_x=params.x;
        if (field=="X" or field=="Y" or field=="Z")
            for (int dim=0;dim<3;dim++)
                params.x[dim]=mrcz::compute_GLL_nodes<T>(params.header.element_shape[dim]);
        printf("params.method= %s\n",params.method.c_str());
        if (params.method=="mrcz")
            if ((params.header.element_shape[0]==params.header.element_shape[1]) and (params.header.element_shape[1]==params.header.element_shape[2]))
                params.compress_regular(field);
            else
                params.compress_irregular(field);
        else if (params.method=="sz2-0")
        {
            vector<vector<T>> backup_x=params.x;
            for (int dim=0;dim<3;dim++)
                if (params.x[dim].size()>1)
                    for (uint32_t i=0;i<params.x[dim].size();i++)
                        params.x[dim][i]=(T)i/(params.x[dim].size()-1)*2-1;
                else
                    params.x[dim][0]=0;
            sz2_1_compress(params,field);
            params.x=backup_x;
        }
        else if (params.method=="sz2-1")
            sz2_1_compress(params,field);
        else if (params.method=="sz3-1")
        {
            vector<vector<T>> backup_x=params.x;
            for (int dim=0;dim<3;dim++)
                if (params.x[dim].size()>1)
                    for (uint32_t i=0;i<params.x[dim].size();i++)
                        params.x[dim][i]=(T)i/(params.x[dim].size()-1)*2-1;
                else
                    params.x[dim][0]=0;
            sz3_2_compress(params,field);
            params.x=backup_x;
        }
        else if (params.method=="sz3-hohs")
            sz3_2_compress(params,field);
        else if (params.method=="zfp")
            zfp_compress(params,field);
        else if (params.method=="sperr")
            sperr_compress(params,field);
        else
            throw runtime_error("未知的压缩方法: "+params.method);
        if (field=="X" or field=="Y" or field=="Z")
            params.x=backup_x;
    }
        
    double br=0,cr=0;
    string compressed_path=params.data_path+"_"+field+"."+params.method;
    ofstream compressed_file(compressed_path,ios::binary);
    uint8_t word_size_during_compression=sizeof(T);
    compressed_file.write(reinterpret_cast<char*>(&word_size_during_compression),1);
    if (params.data_range_min==params.data_range_max)
    {
        unsigned char is_constant=1;
        compressed_file.write(reinterpret_cast<char*>(&is_constant),1);
        compressed_file.write(reinterpret_cast<char*>(&params.data_range_min),sizeof(T));
        br+=1+sizeof(T);
    }
    else
    {
        unsigned char is_constant=0;
        compressed_file.write(reinterpret_cast<char*>(&is_constant),1);
        uint64_t partition_num=params.partitions.size();
        compressed_file.write(reinterpret_cast<char*>(&partition_num),8);
        compressed_file.write(reinterpret_cast<char*>(&params.data_range_min[field]),sizeof(T));
        compressed_file.write(reinterpret_cast<char*>(&params.data_range_max[field]),sizeof(T));
        br+=9+2*sizeof(T);
        for (auto &partition:params.partitions)
        {
            uint64_t zstd_len=partition.zstd_bs.size();
            compressed_file.write(reinterpret_cast<char*>(&zstd_len),8);
        }
        br+=8*params.partitions.size();
        for (auto &partition:params.partitions)
        {
            compressed_file.write(reinterpret_cast<char*>(partition.zstd_bs.data()),partition.zstd_bs.size());
            br+=partition.zstd_bs.size();
            partition.zstd_bs.clear();
            partition.zstd_bs.shrink_to_fit();
        }
    }
    compressed_file.close();
    br*=((double)8/(params.element_num_to_process*params.header.element_size));
    cr=(8*params.header.word_size)/br;
    printf("\033[91mBR= %f bits/value, CR= %f\n\033[0m",br,cr);
}
