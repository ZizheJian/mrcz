#pragma once

namespace mrcz
{
    template<typename T>
    void param_class<T>::read_header(string task)
    {
        if (format=="nekrs")
            read_header_nekrs(task);
        else if (format=="mfem")
            read_header_mfem(task);
        else if (format=="nektar++")
            read_header_nektarpp(task);
    }

    template<typename T>
    void param_class<T>::read_header_nekrs(string task)
    {
        string header_path;
        if (task=="compress")
            header_path=data_path;
        else
            header_path=data_path+"_header.mrcz";
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
        string term;
        int term_index=0;
        while (header_stream>>term)
        {
            if (term_index==1)
                header.word_size=stoi(term);
            else if (term_index==2)
                header.element_shape.push_back(stoi(term));
            else if (term_index==3)
                header.element_shape.push_back(stoi(term));
            else if (term_index==4)
                header.element_shape.push_back(stoi(term));
            else if (term_index==6)
                header.file_element_num=stoull(term);
            else if (term_index==11)
                header.variables=term;
            term_index++;
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
                field_begin_address["Y"]=last_field_end_address;
                field_begin_address["Z"]=last_field_end_address;
                last_field_end_address+=header.file_element_num*header.element_size*header.word_size*3;
                field_end_address["X"]=last_field_end_address;
                field_end_address["Y"]=last_field_end_address;
                field_end_address["Z"]=last_field_end_address;
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
        {
            element_num_to_process=element_end_index-element_begin_index;
            //////// 根据实际需要处理的element范围，调整fields的起始和结束地址 ////////
            for (auto field:fields_to_process)
                if ((field=="X") || (field=="Y") || (field=="Z") || (field=="UX") || (field=="UY") || (field=="UZ"))
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

    template<typename T>
    void param_class<T>::read_header_mfem(string task)
    {
        //////// 读取gf文件的开头 ////////
        string header_path;
        if (task=="compress")
            header_path=data_path;
        else
            header_path=data_path+"_header.mrcz";
        ifstream gf_file(header_path);
        header.original_header="";
        string line;
        while (getline(gf_file,line))
        {
            bool is_empty=true;
            for (char c:line)
                if (!isspace(c))
                {
                    is_empty=false;
                    break;
                }
            if (is_empty)
                break;
            else
                header.original_header+=line+"\n";
        }
        gf_file.close();
        if (task=="decompress")
        {
            uint32_t lastn=header.original_header.rfind('\n');
            uint32_t last2n=header.original_header.rfind('\n',lastn-1);
            string last_line=header.original_header.substr(last2n+1,lastn-last2n-1);
            if (last_line.rfind("DofNum=",0)==0)
                header.dof_num=stoll(last_line.substr(7));
            header.original_header=header.original_header.substr(0,last2n+1);
            ofstream dummy_decompressed_file(data_path+".mrcz.bin");
            dummy_decompressed_file.write(header.original_header.c_str(),header.original_header.size());
            dummy_decompressed_file.put('\n');
            for (uint64_t i=0;i<header.dof_num;i++)
            {
                string s=to_string(i)+"\n";
                dummy_decompressed_file.write(s.c_str(),s.size());
            }
        }
        //////// 读取mfem的一些信息 ////////
        ifstream mesh_f(mesh_path);
        mfem::Mesh mesh(mesh_f,1,1);
        ifstream solution_f(data_path);
        mfem::GridFunction gf(&mesh,solution_f);
        mfem::FiniteElementSpace *fes=gf.FESpace();
        const mfem::FiniteElement *fe=fes->GetFE(0);
        header.word_size=sizeof(mfem::real_t);
        printf("header.word_size=sizeof(mfem::real_t)=%d\n",header.word_size);
        uint32_t order=fe->GetOrder();
        int geomerty=mesh.GetElementBaseGeometry(0);
        if (geomerty==mfem::Geometry::SEGMENT)
            header.element_shape={order+1,1,1};
        else if (geomerty==mfem::Geometry::SQUARE)
            header.element_shape={order+1,order+1,1};
        else if (geomerty==mfem::Geometry::CUBE)
            header.element_shape={order+1,order+1,order+1};
        else
            throw runtime_error("不支持的element形状");
        printf("header.element_shape: %d %d %d\n",header.element_shape[0],header.element_shape[1],header.element_shape[2]);
        header.element_size=fe->GetDof();
        printf("header.element_size=fe->GetDof()=%d\n",header.element_size);
        if (header.element_size!=header.element_shape[0]*header.element_shape[1]*header.element_shape[2])
            throw runtime_error("element_size与element_shape不匹配");
        header.file_element_num=mesh.GetNE();
        printf("header.file_element_num=mesh.GetNE()=%ld\n",header.file_element_num);
        header.variables="";
        header.passive_scalar_num=0;
        header.endian_mode=0;
        header.dof_num=fes->GetNDofs();
        printf("header.dof_num=fes->GetNDofs()=%ld\n",header.dof_num);
        mesh_f.close();
        solution_f.close();
        //////// 生成all_fields_in_file,fields_to_process ////////
        all_fields_in_file.push_back("no_name_field");
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
    }

    template<typename T>
    void param_class<T>::read_header_nektarpp(string task)
    {
        string header_path;
        if (task=="compress")
            header_path=data_path;
        else
            header_path=data_path+"_header.mrcz";
        //////// 读取header文件 ////////
        header.header_doc=new tinyxml2::XMLDocument();
        header.header_doc->LoadFile(header_path.c_str());
        auto *UnstructuredGrid_node=header.header_doc->FirstChildElement("VTKFile")->FirstChildElement("UnstructuredGrid");
        auto *Piece_node=UnstructuredGrid_node->FirstChildElement("Piece");
        uint32_t num_cells=0;
        Piece_node->QueryUnsignedAttribute("NumberOfPoints",&header.element_size);
        Piece_node->QueryUnsignedAttribute("NumberOfCells",&num_cells);
        double a_dim2=pow(static_cast<double>(header.element_size),1.0/2);
        double a_dim3=pow(static_cast<double>(header.element_size),1.0/3);
        double ratio_dim2=abs(log((a_dim2-1)*(a_dim2-1)/num_cells));
        double ratio_dim3=abs(log((a_dim3-1)*(a_dim3-1)*(a_dim3-1)/num_cells));
        if (ratio_dim2<ratio_dim3)
            header.element_shape={1,static_cast<uint32_t>(round(a_dim2)),static_cast<uint32_t>(round(a_dim2))};
        else
            header.element_shape={static_cast<uint32_t>(round(a_dim3)),static_cast<uint32_t>(round(a_dim3)),static_cast<uint32_t>(round(a_dim3))};
        auto *PointData_node=Piece_node->FirstChildElement("PointData");
        for (auto *DataArray_node=PointData_node->FirstChildElement("DataArray");DataArray_node;DataArray_node=DataArray_node->NextSiblingElement("DataArray"))
        {
            string data_type=DataArray_node->Attribute("type");
            header.word_size=0;
            if ((data_type=="Float32") && (header.word_size<4))
                header.word_size=4;
            else if ((data_type=="Float64") && (header.word_size<8))
                header.word_size=8;
            string data_name=DataArray_node->Attribute("Name");
            if (header.variables.find(data_name)==string::npos)
                header.variables+=data_name;
        }
        header.file_element_num=0;
        for (Piece_node=UnstructuredGrid_node->FirstChildElement("Piece");Piece_node;Piece_node=Piece_node->NextSiblingElement("Piece"))
        {
            PointData_node=Piece_node->FirstChildElement("PointData");
            for (auto *DataArray_node=PointData_node->FirstChildElement("DataArray");DataArray_node;DataArray_node=DataArray_node->NextSiblingElement("DataArray"))
                DataArray_node->DeleteChildren();
            header.file_element_num++;
        }
        //////// 生成all_fields_in_file,fields_to_process ////////
        for (uint32_t i=0;i<header.variables.size();i++)
        {
            char c=header.variables[i];
            all_fields_in_file.push_back(string(1,c));
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
    }
}
