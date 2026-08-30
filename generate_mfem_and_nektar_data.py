import os,tarfile,urllib.request
from third_party.get_max_cpus.py_codes.get_max_cpus import get_max_cpus
from third_party.run_command.run_command import run_command

silent_build=False
build_max_cpus=get_max_cpus(8)

def outputs_are_newer_than_script(output_paths):
    script_mtime=os.path.getmtime(__file__)
    return all(os.path.exists(path) and os.path.getmtime(path)>script_mtime for path in output_paths)

mfem_source_path="third_party/mfem"
mfem_build_path="build/mfem_build"
mfem_install_path="install/mfem_install"
os.makedirs(mfem_build_path,exist_ok=True)
run_command(f"cmake -S \"{mfem_source_path}\" -B \"{mfem_build_path}\" "
            f"-DCMAKE_INSTALL_PREFIX=\"{mfem_install_path}\" -DCMAKE_BUILD_TYPE=Release "
            "-DMFEM_ENABLE_TESTING=OFF -DMFEM_ENABLE_EXAMPLES=OFF "
            "-DMFEM_ENABLE_MINIAPPS=OFF",
            silent=silent_build)
run_command(f"cmake --build \"{mfem_build_path}\" --target ex0 ex1 --parallel {build_max_cpus}",
            silent=silent_build)
run_command(f"cmake --install \"{mfem_build_path}\"",silent=silent_build)

mfem_examples_path=os.path.join(mfem_build_path,"examples")
ex0_output_path=os.path.join(mfem_examples_path,"ex0examples")
ex0_outputs=[os.path.join(ex0_output_path,"sol.gf"),os.path.join(ex0_output_path,"mesh.mesh")]
if not outputs_are_newer_than_script(ex0_outputs):
    os.makedirs(ex0_output_path,exist_ok=True)
    run_command("./ex0 --order 20",cwd=mfem_examples_path)
    os.replace(os.path.join(mfem_examples_path,"sol.gf"),ex0_outputs[0])
    os.replace(os.path.join(mfem_examples_path,"mesh.mesh"),ex0_outputs[1])
ex1_output_path=os.path.join(mfem_examples_path,"ex1examples")
ex1_outputs=[os.path.join(ex1_output_path,"sol.gf"),os.path.join(ex1_output_path,"refined.mesh")]
if not outputs_are_newer_than_script(ex1_outputs):
    os.makedirs(ex1_output_path,exist_ok=True)
    run_command("./ex1 --order 4",cwd=mfem_examples_path)
    os.replace(os.path.join(mfem_examples_path,"sol.gf"),ex1_outputs[0])
    os.replace(os.path.join(mfem_examples_path,"refined.mesh"),ex1_outputs[1])

flow_stability_url="http://doc.nektar.info/tutorials/5.10.0/flow-stability/cylinder/flow-stability-cylinder.tar.gz"
flow_stability_path="third_party/flow-stability-cylinder"
flow_stability_archive_path="third_party/flow-stability-cylinder.tar.gz"
if not os.path.exists(flow_stability_path):
    os.makedirs("third_party",exist_ok=True)
    if not os.path.exists(flow_stability_archive_path):
        urllib.request.urlretrieve(flow_stability_url,flow_stability_archive_path)
    with tarfile.open(flow_stability_archive_path,"r:gz") as archive:
        archive.extractall("third_party")

nektar_source_path="third_party/nektar++"
nektar_build_path="build/nektar++_build"
nektar_install_path="install/nektar++_install"
os.makedirs(nektar_build_path,exist_ok=True)
run_command(f"cmake -S \"{nektar_source_path}\" -B \"{nektar_build_path}\" "
            f"-DCMAKE_INSTALL_PREFIX=\"{nektar_install_path}\" -DCMAKE_BUILD_TYPE=Release "
            "-DNEKTAR_BUILD_DEMOS=OFF -DNEKTAR_BUILD_SOLVERS=OFF "
            "-DNEKTAR_BUILD_UNIT_TESTS=OFF -DNEKTAR_BUILD_TESTS=OFF "
            "-DNEKTAR_BUILD_UTILITIES=ON -DNEKTAR_UTILITY_NEKMESH=OFF "
            "-DNEKTAR_UTILITY_FIELDCONVERT=ON -DNEKTAR_USE_SCOTCH=OFF "
            "-DNEKTAR_USE_METIS=OFF",
            silent=silent_build)
run_command(f"cmake --build \"{nektar_build_path}\" --target FieldConvert --parallel {build_max_cpus}",
            silent=silent_build)
run_command(f"cmake --install \"{nektar_build_path}\"",silent=silent_build)

field_convert_path=os.path.join(nektar_install_path,"bin","FieldConvert")
cylinder_base_path=os.path.join(flow_stability_path,"tutorial","base","Cylinder-Base")
cylinder_base_xml_path=cylinder_base_path+".xml"
cylinder_base_fld_path=cylinder_base_path+".fld"
cylinder_base_vtu_path=cylinder_base_path+".vtu"
cylinder_base_inputs=[cylinder_base_xml_path,cylinder_base_fld_path]
vtu_is_current=(outputs_are_newer_than_script([cylinder_base_vtu_path]) and
                os.path.getmtime(cylinder_base_vtu_path)>max(os.path.getmtime(path) for path in cylinder_base_inputs))
if not vtu_is_current:
    run_command(f"\"{field_convert_path}\" --force-output \"{cylinder_base_xml_path}\" "
                f"\"{cylinder_base_fld_path}\" \"{cylinder_base_vtu_path}\"")
