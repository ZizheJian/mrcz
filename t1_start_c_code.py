import argparse,csv,glob,json,os,shlex
from third_party.get_max_cpus.py_codes.get_max_cpus import get_max_cpus
from third_party.run_command.run_command import run_command
from third_party.csv_tools.py_codes.write_csv_item import write_csv_item

project_root=os.path.dirname(os.path.abspath(__file__))
config_directory=os.path.join(project_root,"configs")

def load_json_tasks(json_path):
    with open(json_path,encoding="utf-8") as config_file:
        tasks=json.load(config_file)
    if not isinstance(tasks,list):
        raise TypeError("JSON顶层必须是task list")
    if not tasks:
        raise ValueError("JSON中没有task")
    if any(not isinstance(task,dict) for task in tasks):
        raise TypeError("JSON中的每个task必须是object")
    return tasks

def resolve_json_path(path_or_name):
    path_or_name=os.path.expanduser(path_or_name)
    if os.path.isabs(path_or_name) or os.path.dirname(path_or_name):
        json_path=os.path.abspath(path_or_name)
    else:
        config_filename=path_or_name if path_or_name.endswith(".json") else path_or_name+".json"
        json_path=os.path.join(config_directory,config_filename)
    if not os.path.isfile(json_path):
        raise FileNotFoundError(f"找不到JSON配置: {json_path}")
    return json_path

def remove_file_if_exists(path):
    if os.path.isfile(path):
        os.remove(path)

def cleanup_test_outputs(data_path,method,selected_fields=None):
    if selected_fields is None:
        compressed_paths=glob.glob(f"{glob.escape(data_path)}_*.{glob.escape(method)}")
    else:
        compressed_paths=[f"{data_path}_{field}.{method}" for field in selected_fields]
    for path in compressed_paths:
        remove_file_if_exists(path)
    remove_file_if_exists(f"{data_path}_header.mrcz")
    remove_file_if_exists(f"{data_path}_element_map.mrcz")
    remove_file_if_exists(f"{data_path}.{method}.bin")

def remove_ansi_color(text):
    return text.replace("\033[91m","").replace("\033[0m","")

def parse_rate_distortion(output):
    results={}
    for line in output.splitlines():
        line=remove_ansi_color(line).strip()
        if not line.startswith("Field= "):
            continue
        segments=line.split(",")
        if len(segments)!=8:
            continue
        field=segments[0].split("=",1)[1].strip()
        br=float(segments[3].split("=",1)[1].strip())
        psnr=float(segments[7].split("=",1)[1].strip().split()[0])
        results[field]=(br,psnr)
    return results

def parse_field_throughput(output,task):
    throughput={}
    marker=f"{task} bandwidth per cpu= "
    for line in output.splitlines():
        line=remove_ansi_color(line).strip()
        if not line.startswith("Field= "):
            continue
        segments=line.split(",")
        if len(segments)!=2 or not segments[1].strip().startswith(marker):
            continue
        field=segments[0].split("=",1)[1].strip()
        throughput[field]=float(segments[1].strip()[len(marker):].split()[0])
    return throughput

def parse_mrcz_stream_bitrates(output):
    stream_bitrates={}
    prefix="MRCZ stream bit rates: Field= "
    for line in output.splitlines():
        line=remove_ansi_color(line).strip()
        if not line.startswith(prefix):
            continue
        segments=line[len(prefix):].split(",")
        if len(segments)!=4:
            continue
        field=segments[0].strip()
        coefficient=float(segments[1].split("=",1)[1].strip())
        residual=float(segments[2].split("=",1)[1].strip())
        total=float(segments[3].split("=",1)[1].strip().split()[0])
        stream_bitrates[field]=(coefficient,residual,total)
    return stream_bitrates

def parse_fdwz_prediction_locality(output):
    results={}
    pivot_stride=None
    prefix="FDWZ prediction locality: "
    for line in output.splitlines():
        line=remove_ansi_color(line).strip()
        if not line.startswith(prefix):continue
        payload=line[len(prefix):]
        if payload.startswith("element-index pivot stride="):
            pivot_stride=int(payload.split("=",1)[1])
            continue
        segments=payload.split(",")
        if len(segments)!=5 or not segments[0].startswith("method="):continue
        method=segments[0].split("=",1)[1]
        results[method]={
            "nrmse":float(segments[1].split("=",1)[1]),
            "nmae":float(segments[2].split("=",1)[1]),
            "entropy":float(segments[3].split("=",1)[1].split()[0]),
            "value_num":int(segments[4].split("=",1)[1]),
        }
    return pivot_stride,results

def write_fdwz_prediction_locality(path,output,field,partition_num,selection_rel_eb):
    pivot_stride,results=parse_fdwz_prediction_locality(output)
    expected_methods=("element-local","element-index")
    missing=[method for method in expected_methods if method not in results]
    if pivot_stride is None or missing:
        raise RuntimeError(f"未能解析FDWZ prediction locality结果: pivot_stride={pivot_stride},missing={missing}")
    os.makedirs(os.path.dirname(path) or ".",exist_ok=True)
    field_names=("method","field","partition_num","pivot_stride","selection_rel_eb","value_num","nrmse","nmae","entropy")
    with open(path,"w",newline="") as file_obj:
        writer=csv.DictWriter(file_obj,fieldnames=field_names,lineterminator="\n")
        writer.writeheader()
        for method in expected_methods:
            writer.writerow({"method":method,"field":field,"partition_num":partition_num,"pivot_stride":pivot_stride,"selection_rel_eb":selection_rel_eb,**results[method]})

def validate_parsed_fields(results,expected_fields,result_name,optional_fields=None):
    optional_fields=set() if optional_fields is None else set(optional_fields)
    missing_fields=[field for field in expected_fields if field not in results and field not in optional_fields]
    if missing_fields:
        raise RuntimeError(f"未能解析{result_name}: {','.join(missing_fields)}")

def add_extra_parameters(command_list,extra_parameters):
    for name,value in extra_parameters.items():
        option=name if name.startswith("-") else "-"+name
        if isinstance(value,bool):
            value=str(value).lower()
        command_list=[f"{command} {shlex.quote(option)} {shlex.quote(str(value))}" for command in command_list]
    return command_list

def prepare_tasks(json_path,task_parameters,silent_build,write_to_csv,cleanup_after_each_test):
    required_keys=["name","data_path","format","field_list","eb_mode","eb_list","method_list","max_cpus_list"]
    if json_path is None:
        missing_keys=[key for key in required_keys if task_parameters[key] is None]
        if missing_keys:
            raise ValueError(f"没有提供JSON配置，且缺少参数: {','.join(missing_keys)}")
        tasks=[{key:value for key,value in task_parameters.items() if value is not None}]
    else:
        if not os.path.isabs(json_path):
            raise ValueError("json_path必须是绝对路径")
        if not os.path.isfile(json_path):
            raise FileNotFoundError(f"找不到JSON配置: {json_path}")
        provided_task_parameters={}
        for key,value in task_parameters.items():
            if value is not None:
                provided_task_parameters[key]=value
        tasks=load_json_tasks(json_path)
        if len(tasks)>1 and provided_task_parameters:
            raise ValueError("多任务JSON不能使用task parameters覆盖")
        if len(tasks)==1:
            tasks[0].update(provided_task_parameters)
    runtime_parameters={"silent_build":silent_build,"write_to_csv":write_to_csv,"cleanup_after_each_test":cleanup_after_each_test}
    prepared_tasks=[]
    for config in tasks:
        missing_keys=[key for key in required_keys if key not in config]
        if missing_keys:
            raise ValueError(f"task缺少参数: {','.join(missing_keys)}")
        if config["format"]=="mfem" and not config.get("mesh_path"):
            raise ValueError(f"{config['name']}使用MFEM格式但没有mesh_path")
        config={**config,**runtime_parameters}
        if not os.path.isabs(config["data_path"]):
            config["data_path"]=os.path.join(project_root,config["data_path"])
        if config.get("mesh_path") is not None and not os.path.isabs(config["mesh_path"]):
            config["mesh_path"]=os.path.join(project_root,config["mesh_path"])
        if config.get("csv_path") is not None and not os.path.isabs(config["csv_path"]):
            config["csv_path"]=os.path.join(project_root,config["csv_path"])
        if config.get("prediction_locality_csv_path") is not None and not os.path.isabs(config["prediction_locality_csv_path"]):
            config["prediction_locality_csv_path"]=os.path.join(project_root,config["prediction_locality_csv_path"])
        prepared_tasks.append(config)
    return prepared_tasks

def build_project(build_max_cpus,silent_build):
    # SPERR的CMakeLists无法由主项目直接安装，需要单独构建和安装。
    sperr_source_path=os.path.join(project_root,"third_party/sperr")
    sperr_build_path=os.path.join(project_root,"build/sperr_build")
    sperr_install_path=os.path.join(project_root,"install/sperr_install")
    os.makedirs(sperr_build_path,exist_ok=True)
    if not os.path.exists(os.path.join(sperr_build_path,"CMakeCache.txt")):
        run_command(f"cmake -S \"{sperr_source_path}\" -B \"{sperr_build_path}\" -DCMAKE_INSTALL_PREFIX=\"{sperr_install_path}\" -DBUILD_SHARED_LIBS=OFF",silent=silent_build)
    run_command(f"cmake --build {sperr_build_path} --parallel {build_max_cpus}",silent=silent_build)
    run_command(f"cmake --install {sperr_build_path}",silent=silent_build)
    mfem_source_path=os.path.join(project_root,"third_party/mfem")
    mfem_build_path=os.path.join(project_root,"build/mfem_build")
    mfem_install_path=os.path.join(project_root,"install/mfem_install")
    os.makedirs(mfem_build_path,exist_ok=True)
    if not os.path.exists(os.path.join(mfem_build_path,"CMakeCache.txt")):
        run_command(f"cmake -S \"{mfem_source_path}\" -B \"{mfem_build_path}\" -DCMAKE_INSTALL_PREFIX=\"{mfem_install_path}\"",silent=silent_build)
    run_command(f"cmake --build {mfem_build_path} --parallel {build_max_cpus}",silent=silent_build)
    run_command(f"cmake --install {mfem_build_path}",silent=silent_build)
    main_source_path=os.path.join(project_root,"c_codes")
    main_build_path=os.path.join(project_root,"build")
    main_install_path=os.path.join(project_root,"install")
    os.makedirs(main_build_path,exist_ok=True)
    run_command(f"cmake -S \"{main_source_path}\" -B \"{main_build_path}\" -DCMAKE_INSTALL_PREFIX=\"{main_install_path}\"",silent=silent_build)
    run_command(f"cmake --build {main_build_path} --parallel {build_max_cpus}",silent=silent_build)
    run_command(f"cmake --install {main_build_path}",silent=silent_build)

def run_experiment(json_path=None,name=None,data_path=None,mesh_path=None,format=None,field_list=None,eb_mode=None,eb_list=None,
                   method_list=None,max_cpus_list=None,partition_num=None,elements=None,extra_parameters=None,csv_path=None,
                   build_before_test=True,build_max_cpus=None,silent_build=False,write_to_csv=True,cleanup_after_each_test=False):
    task_parameters={"name":name,"data_path":data_path,"mesh_path":mesh_path,"format":format,"field_list":field_list,
                     "eb_mode":eb_mode,"eb_list":eb_list,"method_list":method_list,"max_cpus_list":max_cpus_list,
                     "partition_num":partition_num,"elements":elements,"extra_parameters":extra_parameters,"csv_path":csv_path}
    prepared_tasks=prepare_tasks(json_path,task_parameters,silent_build,write_to_csv,cleanup_after_each_test)
    if build_before_test:
        build_project(get_max_cpus(8) if build_max_cpus is None else build_max_cpus,silent_build)
    for config in prepared_tasks:
        task_name=config["name"]
        write_to_csv=config.get("write_to_csv",True)
        cleanup_after_each_test=config.get("cleanup_after_each_test",False)
        format=config["format"]
        data_path=config["data_path"]
        mesh_path=config.get("mesh_path")
        fields=config["field_list"]
        eb_mode=config["eb_mode"]
        eb_list=config["eb_list"]
        method_list=config["method_list"]
        max_cpus_list=config["max_cpus_list"]
        partition_num=config.get("partition_num")
        elements=config.get("elements")
        extra_parameters=config.get("extra_parameters") or {}
        csv_path=config.get("csv_path")
        prediction_locality_csv_path=config.get("prediction_locality_csv_path")
        compression_only=config.get("compression_only",False)
        trials=config.get("trials",1)
        if not isinstance(trials,int) or trials<=0:
            raise ValueError(f"{task_name}的trials必须是正整数")
        csv_file_prefix=data_path
        if csv_path is not None:
            if write_to_csv:
                os.makedirs(csv_path,exist_ok=True)
            csv_file_prefix=os.path.join(csv_path,os.path.basename(data_path))

        for method in method_list:
            for eb in eb_list:
                for max_cpus in max_cpus_list:
                    for trial in range(1,trials+1):
                        command_list=[os.path.join(project_root,"install/bin/mrcz_compress"),os.path.join(project_root,"install/bin/mrcz_decompress")]
                        if format in ["nekrs","nektar++"]:
                            command_list=[f"{shlex.quote(command)} -format {shlex.quote(format)} -i {shlex.quote(data_path)}" for command in command_list]
                        elif format=="mfem":
                            command_list=[f"{shlex.quote(command)} -format {shlex.quote(format)} -i {shlex.quote(data_path)} -mesh {shlex.quote(mesh_path)}" for command in command_list] # type: ignore
                        command_list=[f"{command} -M {shlex.quote(eb_mode)} {eb} -method {shlex.quote(method)}" for command in command_list]
                        command_list=[f"{command} -max-cpus {max_cpus}" for command in command_list]
                        if partition_num is not None:
                            command_list=[f"{command} -partition-num {partition_num}" for command in command_list]
                        if elements is not None:
                            command_list=[f"{command} -elements {shlex.quote(str(elements))}" for command in command_list]
                        if fields is not None:
                            fields_argument=shlex.quote("["+",".join(fields)+"]")
                            command_list=[f"{command} -fields {fields_argument}" for command in command_list] # type: ignore
                        command_list=add_extra_parameters(command_list,extra_parameters)
                        compression_output=run_command(command_list[0])
                        if prediction_locality_csv_path is not None:
                            if len(fields)!=1 or eb_mode!="REL":
                                raise ValueError("prediction locality CSV要求单个field和REL EB")
                            write_fdwz_prediction_locality(prediction_locality_csv_path,compression_output,fields[0],partition_num,eb)
                        if compression_only:
                            continue
                        decompression_output=run_command(command_list[1])
                        if write_to_csv:
                            compression_throughput=parse_field_throughput(compression_output,"compression")
                            decompression_throughput=parse_field_throughput(decompression_output,"decompression")
                            stream_bitrates=parse_mrcz_stream_bitrates(compression_output)
                            rate_distortion=parse_rate_distortion(decompression_output)
                            constant_fields={field for field in fields if f"Field {field} is compressed as constant:" in compression_output}
                            result_context=f"{task_name}/{method}/EB={eb}/max_cpus={max_cpus}/trial={trial}"
                            validate_parsed_fields(compression_throughput,fields,f"{result_context}的compression throughput",constant_fields)
                            validate_parsed_fields(decompression_throughput,fields,f"{result_context}的decompression throughput",constant_fields)
                            validate_parsed_fields(rate_distortion,fields,f"{result_context}的rate-distortion results")
                            if extra_parameters.get("analyze-mrcz-stream-sizes"):
                                validate_parsed_fields(stream_bitrates,fields,f"{result_context}的MRCZ stream bit rates")
                            csv_method_name=method
                            for field,(br,psnr) in rate_distortion.items():
                                result_row={f"{eb_mode} EB":eb}
                                if trials>1:
                                    result_row["Trial"]=trial
                                rate_distortion_path=f"{csv_file_prefix}_{field}.csv"
                                write_csv_item(rate_distortion_path,result_row,f"{csv_method_name}_BR",br)
                                write_csv_item(rate_distortion_path,result_row,f"{csv_method_name}_PSNR",psnr)
                                if field in stream_bitrates:
                                    coefficient_br,residual_br,separate_total_br=stream_bitrates[field]
                                    write_csv_item(rate_distortion_path,result_row,f"{csv_method_name}_Coefficient_BR",coefficient_br)
                                    write_csv_item(rate_distortion_path,result_row,f"{csv_method_name}_Residual_BR",residual_br)
                                    write_csv_item(rate_distortion_path,result_row,f"{csv_method_name}_Separately_Compressed_BR",separate_total_br)

                                throughput_path=f"{csv_file_prefix}_{field}_throughput.csv"
                                throughput_row={**result_row,"Max CPUs":max_cpus}
                                if partition_num is not None:
                                    throughput_row["Partition Num"]=partition_num
                                if field in compression_throughput:
                                    write_csv_item(throughput_path,throughput_row,
                                        f"{csv_method_name} Compression Throughput",compression_throughput[field])
                                if field in decompression_throughput:
                                    write_csv_item(throughput_path,throughput_row,
                                        f"{csv_method_name} Decompression Throughput",decompression_throughput[field])
                        if cleanup_after_each_test:
                            cleanup_test_outputs(data_path,method,fields)

if __name__=="__main__":
    parser=argparse.ArgumentParser(description="运行configs目录中的MRCZ实验配置")
    parser.add_argument("json",help="configs中的JSON配置名称，或JSON文件路径")
    arguments=parser.parse_args()
    try:
        json_path=resolve_json_path(arguments.json)
    except FileNotFoundError as error:
        parser.error(str(error))
    run_experiment(json_path=json_path)
