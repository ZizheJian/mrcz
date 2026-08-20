#pragma once

#include "mrcz.hpp"

template<typename T>
void read_data_nekrs(mrcz::param_class<T> &params,string field);
template<typename T>
void read_data_mfem(mrcz::param_class<T> &params,string field);
template<typename T>
void read_data_nektarpp(mrcz::param_class<T> &params,string field);
template<typename T>
void write_data_nekrs(mrcz::param_class<T> &params,string field);
template<typename T>
void write_data_mfem(mrcz::param_class<T> &params,string field);
template<typename T>
void write_data_nektarpp(mrcz::param_class<T> &params,string field);

template<typename T>
void read_data(mrcz::param_class<T> &params,string field)//分配好params.data的空间后，调用这个函数来读取指定field的数据到params.data中
{
    if (params.format=="nekrs")
        read_data_nekrs(params,field);
    else if (params.format=="mfem")
        read_data_mfem(params,field);
    else if (params.format=="nektar++")
        read_data_nektarpp(params,field);
    else
        throw runtime_error("不支持的格式");
}

template<typename T>
void read_data_nekrs(mrcz::param_class<T> &params,string field)
{
    auto &header=params.header;
    ifstream data_file(params.data_path,ios::binary);
    data_file.seekg(params.field_begin_address[field],ios::beg);
    if (field=="X" or field=="Y" or field=="Z" or field=="UX" or field=="UY" or field=="UZ")
    {
        vector<T> temp_data(params.element_num_to_process*header.element_size*3);
        data_file.read(reinterpret_cast<char*>(temp_data.data()),params.element_num_to_process*header.element_size*3*header.word_size);
        if (field=="X" or field=="UX")
            for (uint64_t element_id=0;element_id<params.element_num_to_process;element_id++)
                memcpy(&params.data[element_id*header.element_size],&temp_data[3*element_id*header.element_size],
                    header.element_size*header.word_size);
        if (field=="Y" or field=="UY")
            for (uint64_t element_id=0;element_id<params.element_num_to_process;element_id++)
                memcpy(&params.data[element_id*header.element_size],&temp_data[(3*element_id+1)*header.element_size],
                    header.element_size*header.word_size);
        if (field=="Z" or field=="UZ")
            for (uint64_t element_id=0;element_id<params.element_num_to_process;element_id++)
                memcpy(&params.data[element_id*header.element_size],&temp_data[(3*element_id+2)*header.element_size],
                    header.element_size*header.word_size);
    }
    else
        data_file.read(reinterpret_cast<char*>(params.data.data()),params.element_num_to_process*header.element_size*header.word_size);
    data_file.close();
}

template<typename T>
void read_data_mfem(mrcz::param_class<T> &params,string field)
{
    auto &header=params.header;
    ifstream data_file(params.data_path);
    ifstream mesh_file(params.mesh_path);
    mfem::Mesh mesh(mesh_file,1,1);
    mfem::GridFunction gf(&mesh,data_file);
    mfem::FiniteElementSpace *fes=gf.FESpace();
    vector<T> neg_coordinate_quant_width(3);
    vector<unordered_map<int32_t,uint32_t>> coordinate_maps(3);
    for (int dim=0;dim<3;dim++)
    {
        if (params.x[dim].size()>1)
            neg_coordinate_quant_width[dim]=1.1/(params.x[dim][1]-params.x[dim][0])*2;
        else
            neg_coordinate_quant_width[dim]=0;
        for (uint64_t i=0;i<params.x[dim].size();i++)
            coordinate_maps[dim][static_cast<int32_t>(round((1+params.x[dim][i])/2*neg_coordinate_quant_width[dim]))]=i;
    }
    for (uint64_t e=0;e<params.element_num_to_process;e++)
    {
        const mfem::FiniteElement *fe=fes->GetFE(e);
        if (fe->GetOrder()!=fes->GetFE(0)->GetOrder())
            throw runtime_error("不支持不同order的element");
        if (mesh.GetElementBaseGeometry(e)!=mesh.GetElementBaseGeometry(0))
            throw runtime_error("不支持不同形状的element");
        if (fe->GetDof()!=fes->GetFE(0)->GetDof())
            throw runtime_error("不支持不同dof的element");
        const mfem::IntegrationRule &ir=fe->GetNodes();
        mfem::Array<int> vdofs;
        fes->GetElementVDofs(e,vdofs);
        for (uint32_t i=0;i<header.element_size;i++)
        {
            const mfem::IntegrationPoint &ip=ir.IntPoint(i);
            int idx=coordinate_maps[0][static_cast<int32_t>(round(ip.x*neg_coordinate_quant_width[0]))]*params.x[1].size()*params.x[2].size()+
                    coordinate_maps[1][static_cast<int32_t>(round(ip.y*neg_coordinate_quant_width[1]))]*params.x[2].size()+
                    coordinate_maps[2][static_cast<int32_t>(round(ip.z*neg_coordinate_quant_width[2]))];
            params.data[e*header.element_size+idx]=gf(vdofs[i]);
        }
    }
    mesh_file.close();
    data_file.close();
}

template<typename T>
void read_data_nektarpp(mrcz::param_class<T> &params,string field)
{
    tinyxml2::XMLDocument data_doc;
    data_doc.LoadFile(params.data_path.c_str());
    auto *UnstructuredGrid_node=data_doc.FirstChildElement("VTKFile")->FirstChildElement("UnstructuredGrid");
    uint64_t element_id=0;
    for (auto *Piece_node=UnstructuredGrid_node->FirstChildElement("Piece");Piece_node;Piece_node=Piece_node->NextSiblingElement("Piece"))
    {
        auto *PointData_node=Piece_node->FirstChildElement("PointData");
        for (auto *DataArray_node=PointData_node->FirstChildElement("DataArray");DataArray_node;DataArray_node=DataArray_node->NextSiblingElement("DataArray"))
        {
            if (string(DataArray_node->Attribute("Name"))==field)
            {
                char *text=(char*)(DataArray_node->GetText());
                for (uint32_t point_id=0;point_id<params.header.element_size;point_id++)
                {
                    char *endptr;
                    params.data[element_id*params.header.element_size+point_id]=static_cast<T>(strtod(text,&endptr)); 
                    text=endptr;
                }
            }
        }
        element_id++;
    }
}

template<typename T>
void write_data(mrcz::param_class<T> &params,string field)
{
    assert(params.header.word_size==sizeof(T));
    if (params.format=="nekrs")
        write_data_nekrs(params,field);
    else if (params.format=="mfem")
        write_data_mfem(params,field);
    else if (params.format=="nektar++")
        write_data_nektarpp(params,field);
    else
        throw runtime_error("不支持的格式");
}

template<typename T>
void write_data_nekrs(mrcz::param_class<T> &params,string field)
{
    fstream decompressed_file(params.data_path+"."+params.method+".bin",ios::app|ios::binary);
    if (field=="UX" || field=="UY" || field=="UZ")
    {
        if (field=="UX")
        {
            decompressed_file.seekp(0,ios::end);
            // assert(decompressed_file.tellp()==params.field_begin_address[field]);
            vector<T> temp_data(params.element_num_to_process*params.header.element_size*3,0);
            for (uint64_t element_id=0;element_id<params.element_num_to_process;element_id++)
                memcpy(&temp_data[element_id*params.header.element_size*3],&params.decompressed_data[element_id*params.header.element_size],params.header.element_size*sizeof(T));
            decompressed_file.write(reinterpret_cast<char*>(temp_data.data()),params.element_num_to_process*params.header.element_size*3*sizeof(T));
        }
        else if (field=="UY")
        {
            decompressed_file.seekp(0,ios::end);
            // assert(decompressed_file.tellp()==params.field_end_address[field]);
            vector<T> temp_data(params.element_num_to_process*params.header.element_size*3,0);
            decompressed_file.seekg(params.field_begin_address[field]);
            decompressed_file.read(reinterpret_cast<char*>(temp_data.data()),params.element_num_to_process*params.header.element_size*3*sizeof(T));
            for (uint64_t element_id=0;element_id<params.element_num_to_process;element_id++)
                memcpy(&temp_data[(3*element_id+1)*params.header.element_size],&params.decompressed_data[element_id*params.header.element_size],
                    params.header.element_size*sizeof(T));
            decompressed_file.write(reinterpret_cast<char*>(temp_data.data()),params.element_num_to_process*params.header.element_size*3*sizeof(T));
        }
        else if (field=="UZ")
        {
            decompressed_file.seekp(0,ios::end);
            // assert(decompressed_file.tellp()==params.field_end_address[field]);
            vector<T> temp_data(params.element_num_to_process*params.header.element_size*3,0);
            decompressed_file.seekg(params.field_begin_address[field]);
            decompressed_file.read(reinterpret_cast<char*>(temp_data.data()),params.element_num_to_process*params.header.element_size*3*sizeof(T));
            for (uint64_t element_id=0;element_id<params.element_num_to_process;element_id++)
                memcpy(&temp_data[(3*element_id+2)*params.header.element_size],&params.decompressed_data[element_id*params.header.element_size],
                    params.header.element_size*sizeof(T));
            decompressed_file.write(reinterpret_cast<char*>(temp_data.data()),params.element_num_to_process*params.header.element_size*3*sizeof(T));
        }
    }
    else
    {
        decompressed_file.seekp(0,ios::end);
        // assert(decompressed_file.tellp()==params.field_begin_address[field]);
        decompressed_file.write(reinterpret_cast<char*>(params.decompressed_data.data()),params.element_num_to_process*params.header.element_size*params.header.word_size);
    }
    decompressed_file.close();
}

template<typename T>
void write_data_mfem(mrcz::param_class<T> &params,string field)
{
    // 先建立一个dof的值的vector和dof重复次数的vector
    // 然后读取mesh文件和假的解压缩文件，确认每一个解压出的数据对应解压文件里的哪一个dof
    // 然后把解压数据累加到对应的dof上
    // 最后除以重复次数，写出解压缩文件
    auto &header=params.header;
    vector<T> dof_values(header.dof_num,0);
    vector<uint64_t> dof_duplicate_times(header.dof_num,0);
    ifstream data_file(params.data_path+"."+params.method+".bin");
    ifstream mesh_file(params.mesh_path);
    mfem::Mesh mesh(mesh_file,1,1);
    mfem::GridFunction gf(&mesh,data_file);
    mfem::FiniteElementSpace *fes=gf.FESpace();
    vector<T> neg_coordinate_quant_width(3);
    vector<unordered_map<int32_t,uint32_t>> coordinate_maps(3);
    for (int dim=0;dim<3;dim++)
    {
        if (params.x[dim].size()>1)
            neg_coordinate_quant_width[dim]=1.1/(params.x[dim][1]-params.x[dim][0])*2;
        else
            neg_coordinate_quant_width[dim]=0;
        for (uint64_t i=0;i<params.x[dim].size();i++)
            coordinate_maps[dim][static_cast<int32_t>(round((1+params.x[dim][i])/2*neg_coordinate_quant_width[dim]))]=i;
    }
    for (uint64_t e=0;e<params.element_num_to_process;e++)
    {
        const mfem::FiniteElement *fe=fes->GetFE(e);
        const mfem::IntegrationRule &ir=fe->GetNodes();
        mfem::Array<int> vdofs;
        fes->GetElementVDofs(e,vdofs);
        for (uint32_t i=0;i<header.element_size;i++)
        {
            const mfem::IntegrationPoint &ip=ir.IntPoint(i);
            int idx=coordinate_maps[0][static_cast<int32_t>(round(ip.x*neg_coordinate_quant_width[0]))]*params.x[1].size()*params.x[2].size()+
                    coordinate_maps[1][static_cast<int32_t>(round(ip.y*neg_coordinate_quant_width[1]))]*params.x[2].size()+
                    coordinate_maps[2][static_cast<int32_t>(round(ip.z*neg_coordinate_quant_width[2]))];
            // e*header.element_size+idx刚解压出的一个数，gf(vdofs[i])是这个数对应的假的dof的值，同时也是dof编号
            // printf("%f %lu\n",gf(vdofs[i]),uint64_t(gf(vdofs[i])));
            // getchar();
            dof_values[uint64_t(gf(vdofs[i]))]+=params.decompressed_data[e*header.element_size+idx];
            dof_duplicate_times[uint64_t(gf(vdofs[i]))]++;
        }
    }
    mesh_file.close();
    data_file.close();
    for (uint64_t i=0;i<header.dof_num;i++)
    {
        assert(dof_duplicate_times[i]>0);
        dof_values[i]/=dof_duplicate_times[i];
    }
    ofstream decompressed_file(params.data_path+"."+params.method+".bin");
    decompressed_file.write(header.original_header.c_str(),header.original_header.size());
    decompressed_file.put('\n');
    char buf[64];
    for (uint64_t i=0;i<header.dof_num;i++)
    {
        auto [ptr,ec]=to_chars(buf,buf+64,dof_values[i],chars_format::general,16);
        *ptr++='\n';
        decompressed_file.write(buf,ptr-buf);
    }
    decompressed_file.close();
}

template<typename T>
void write_data_nektarpp(mrcz::param_class<T> &params,string field)
{
    // Implementation for writing data in nektar++ format
}
