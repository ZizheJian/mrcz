#include "main.hpp"

template<typename T>
unordered_map<string,vector<double>> decompress1(mrcz::param_class<T> &params);
template<typename T>
vector<double> decompress2(mrcz::param_class<T> &params,string field);

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
    
    if (params.format.empty())
        throw runtime_error("必须指定数据格式，使用-format参数");
    if (params.data_path.empty())
        throw runtime_error("必须指定原始数据路径，使用-i参数。如果原始数据不存在，寻找解压产生的独立的header文件，文件名以_header.mrcz为结尾），去掉这一后缀作为数据路径");
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
    
    params.read_header("decompress");
    params.default_init();
    params.save_header_and_element_map("decompress");
    params.get_max_cpus();
    
    unordered_map<string,vector<double>> result;
    if (params.header.word_size==8)
    {
        mrcz::param_class<double> params_double(params);
        result=decompress1(params_double);
        for (auto field:params.fields_to_process)
            params.abs_eb[field]=params_double.abs_eb[field];
    }
    else
        result=decompress1(params);
}

template<typename T>
unordered_map<string,vector<double>> decompress1(mrcz::param_class<T> &params) //对不同的field调用decompress2，并修改T至实际压缩的类型
{
    unordered_map<string,vector<double>> result;
    for (auto field:params.fields_to_process)
    {
        printf("Field= %s\n",field.c_str());
        if (field=="X")
        {
            auto header=params.header;
            ifstream compressed_file(params.data_path+"_"+field+".mrcz",ios::binary);
            vector<T> coordinate_data(header.file_element_num*3*header.element_size);
            compressed_file.read(reinterpret_cast<char*>(coordinate_data.data()),header.file_element_num*3*header.element_size*header.word_size);
            compressed_file.close();
            ofstream decompressed_file(params.data_path+".mrcz.bin",ios::app|ios::binary);
            decompressed_file.seekp(0,ios::end);
            decompressed_file.write(reinterpret_cast<char*>(coordinate_data.data()),header.file_element_num*3*header.element_size*header.word_size);
            decompressed_file.close();
            continue;
        }
        ifstream compressed_file(params.data_path+"_"+field+".mrcz",ios::binary);
        uint8_t word_size_during_compression;
        compressed_file.read(reinterpret_cast<char*>(&word_size_during_compression),1);
        compressed_file.close();
        if ((is_same<T,float>::value) and (word_size_during_compression==8)) //原始数据类型是float，但是压缩时使用了double
        {
            mrcz::param_class<double> params_double(params);
            result[field]=decompress2(params_double,field);
            params.abs_eb[field]=params_double.abs_eb[field];
        }
        else
        {
            result[field]=decompress2(params,field);
        }
    }
    return result;
}

template<typename T>
vector<double> decompress2(mrcz::param_class<T> &params,string field)
{
    ifstream compressed_file(params.data_path+"_"+field+".mrcz",ios::binary);
    compressed_file.seekg(1,ios::beg);
    uint64_t block_num;
    compressed_file.read(reinterpret_cast<char*>(&block_num),8);
    params.blocks.resize(block_num);
    for (uint64_t i=0;i<block_num;i++)
    {
        auto &block=params.blocks[i];
        block.block_id=i;
        block.begin_index=(i*params.element_num_to_process)/block_num;
        block.end_index=((i+1)*params.element_num_to_process)/block_num;
        block.element_num=block.end_index-block.begin_index;
    }
    params.data_range_min[field]=0;
    params.data_range_max[field]=0;
    if (is_same<T,float>::value)
    {
        compressed_file.read(reinterpret_cast<char*>(&params.data_range_min[field]),sizeof(T));
        compressed_file.read(reinterpret_cast<char*>(&params.data_range_max[field]),sizeof(T));
    }
    else
    {
        float data_range_min_float,data_range_max_float;
        compressed_file.read(reinterpret_cast<char*>(&data_range_min_float),sizeof(float));
        compressed_file.read(reinterpret_cast<char*>(&data_range_max_float),sizeof(float));
        params.data_range_min[field]=data_range_min_float;
        params.data_range_max[field]=data_range_max_float;
    }
    for (uint64_t i=0;i<block_num;i++)
    {
        auto &block=params.blocks[i];
        uint64_t zstd_len;
        compressed_file.read(reinterpret_cast<char*>(&zstd_len),8);
        block.zstd_bs.resize(zstd_len);
    }
    for (uint64_t i=0;i<block_num;i++)
    {
        auto &block=params.blocks[i];
        compressed_file.read(reinterpret_cast<char*>(block.zstd_bs.data()),block.zstd_bs.size());
    }
    double br=((double)compressed_file.tellg()*8)/(params.element_num_to_process*params.header.element_size);
    double cr=(params.header.word_size*8)/br;
    compressed_file.close();

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

    printf("params.method= %s\n",params.method.c_str());
    if (params.method=="mrcz")
        if ((params.header.element_shape[0]==params.header.element_shape[1]) and (params.header.element_shape[1]==params.header.element_shape[2]))
            params.decompress_regular(field);
        else
            params.decompress_irregular(field);
    else if (params.method=="sz2-0")
    {
        vector<vector<T>> backup_x=params.x;
        for (int dim=0;dim<3;dim++)
            if (params.x[dim].size()>1)
                for (uint32_t i=0;i<params.x[dim].size();i++)
                    params.x[dim][i]=(T)i/(params.x[dim].size()-1)*2-1;
            else
                params.x[dim][0]=0;
        sz2_1_decompress(params,field);
        params.x=backup_x;
    }
    else if (params.method=="sz2-1")
        sz2_1_decompress(params,field);
    else if (params.method=="sz3-1")
    {
        vector<vector<T>> backup_x=params.x;
        for (int dim=0;dim<3;dim++)
            if (params.x[dim].size()>1)
                for (uint32_t i=0;i<params.x[dim].size();i++)
                    params.x[dim][i]=(T)i/(params.x[dim].size()-1)*2-1;
            else
                params.x[dim][0]=0;
        sz3_2_decompress(params,field);
        params.x=backup_x;
    }
    else if (params.method=="sz3-hohs")
        sz3_2_decompress(params,field);
    else if (params.method=="zfp")
        zfp_decompress(params,field);
    else if (params.method=="sperr")
        sperr_decompress(params,field);
    else
        throw runtime_error("未知的解压方法: "+params.method);

    //清空zstd_bs
    for (auto &block:params.blocks)
    {
        block.zstd_bs.clear();
        block.zstd_bs.shrink_to_fit();
    }
    
    //写出解压数据
    if (params.fields_to_process.size()!=params.all_fields_in_file.size() || params.header.file_element_num!=params.element_num_to_process)
        printf("没有选择全部的field，或者没有选择全部的element进行解压，不生成解压文件\n");
    else
        if (is_same<T,double>::value and params.header.word_size==4) 
        {
            mrcz::param_class<float> params_float(params);
            params_float.decompressed_data.resize(params.decompressed_data.size());
            #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
            for (uint64_t block_id=0;block_id<params.blocks.size();block_id++)
            {
                auto &block=params.blocks[block_id];
                uint64_t begin_address=block.begin_index*params.header.element_size;
                uint64_t end_address=block.end_index*params.header.element_size;
                for (uint64_t i=begin_address;i<end_address;i++)
                    params_float.decompressed_data[i]=params.decompressed_data[i];
            }
            write_data(params_float,field);
            params_float.decompressed_data.clear();
            params_float.decompressed_data.shrink_to_fit();
        }
        else
            write_data(params,field);
    
    //计算压缩质量
    vector<double> result;
    ifstream data_file(params.data_path,ios::binary);
    if (data_file.is_open())
    {
        data_file.close();
        params.data.resize(params.element_num_to_process*params.header.element_size);
        if ((params.header.word_size==4) and (is_same<T,double>::value))
        {
            mrcz::param_class<float> params_float(params);
            params_float.data.resize(params.element_num_to_process*params.header.element_size);
            read_data(params_float,field);
            #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
            for (uint64_t block_id=0;block_id<params.blocks.size();block_id++)
            {
                auto &block=params.blocks[block_id];
                uint64_t begin_address=block.begin_index*params.header.element_size;
                uint64_t end_address=block.end_index*params.header.element_size;
                for (uint64_t i=begin_address;i<end_address;i++)
                    params.data[i]=params_float.data[i];
            }
            params_float.data.clear();
            params_float.data.shrink_to_fit();
        }
        else
            read_data(params,field);
        double max_err=0,mse=0,psnr=0;
        #pragma omp parallel for reduction(max:max_err) reduction(+:mse) schedule(static) num_threads(params.max_cpus)
        for (uint64_t block_id=0;block_id<params.blocks.size();block_id++)
        {
            auto &block=params.blocks[block_id];
            uint64_t begin_address=block.begin_index*params.header.element_size;
            uint64_t end_address=block.end_index*params.header.element_size;
            for (uint64_t i=begin_address;i<end_address;i++)
            {
                T err=fabs(params.data[i]-params.decompressed_data[i]);
                if (err>params.abs_eb[field]*1.1)
                {
                    printf("index= %lu, original= %f, decompressed= %f, err= %f, abs_eb= %f\n",i,params.data[i],params.decompressed_data[i],err,params.abs_eb[field]);
                    exit(0);
                }
                max_err=max(max_err,(double)err);
                mse+=err*err;
            }
        }
        params.data.clear();
        params.data.shrink_to_fit();
        mse/=(params.element_num_to_process*params.header.element_size);
        psnr=10*log10(pow(params.data_range_max[field]-params.data_range_min[field],2)/mse);
        result={br,cr,max_err,mse,psnr};
    }
    else
    {
        printf("原始数据不存在，无法计算压缩质量\n");
        result={br,cr,0,0,0};
    }
    params.decompressed_data.clear();
    params.decompressed_data.shrink_to_fit();
    printf("\033[91mField= %s, abs_eb= %.15f, BR= %.15f, CR= %.15f, max_err= %.15f (%f%%), MSE= %.15f, PSNR= %.15f dB\033[0m\n",
        field.c_str(),params.abs_eb[field],result[0],result[1],result[2],result[2]/params.abs_eb[field]*100,result[3],result[4]);
    return result;
}