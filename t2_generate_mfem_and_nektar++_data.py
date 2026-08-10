import os,subprocess,tarfile,urllib.request
from third_party.get_max_cpus.py_codes.get_max_cpus import get_max_cpus
from third_party.run_command.run_command import run_command

redbegin="\033[91m"
redend="\033[0m"
silent_build=False
build_max_cpus=get_max_cpus(8)

######## build mfem ########
mfem_source_path="third_party/mfem"
mfem_build_path="build/mfem_build"
mfem_install_path="install/mfem_install"
os.makedirs(mfem_build_path,exist_ok=True)
if not os.path.exists(os.path.join(mfem_build_path,"CMakeCache.txt")):
    run_command(f"cmake -S \"{mfem_source_path}\" -B \"{mfem_build_path}\" -DCMAKE_INSTALL_PREFIX=\"{mfem_install_path}\"",
                silent=silent_build)
run_command(f"cmake --build {mfem_build_path} --parallel {build_max_cpus}",silent=silent_build)
run_command(f"cmake --build {mfem_build_path} --target examples --parallel {build_max_cpus}",silent=silent_build)
run_command(f"cmake --install {mfem_build_path}",silent=silent_build)

######## generate mfem data ########
os.makedirs("build/mfem_build/examples/ex0examples",exist_ok=True)
run_command("./ex0 --order 20",cwd="build/mfem_build/examples")
os.rename("build/mfem_build/examples/sol.gf","build/mfem_build/examples/ex0examples/sol.gf")
os.rename("build/mfem_build/examples/mesh.mesh","build/mfem_build/examples/ex0examples/mesh.mesh")
os.makedirs("build/mfem_build/examples/ex1examples",exist_ok=True)
run_command("./ex1 --order 4",cwd="build/mfem_build/examples")
os.rename("build/mfem_build/examples/sol.gf","build/mfem_build/examples/ex1examples/sol.gf")
os.rename("build/mfem_build/examples/refined.mesh","build/mfem_build/examples/ex1examples/refined.mesh")

######## download flow-stability-cylinder ########
flow_stability_url="http://doc.nektar.info/tutorials/5.10.0/flow-stability/cylinder/flow-stability-cylinder.tar.gz"
flow_stability_path="third_party/flow-stability-cylinder"
flow_stability_archive_path="third_party/flow-stability-cylinder.tar.gz"
if not os.path.exists(flow_stability_path):
    os.makedirs("third_party",exist_ok=True)
    if not os.path.exists(flow_stability_archive_path):
        urllib.request.urlretrieve(flow_stability_url,flow_stability_archive_path)
    with tarfile.open(flow_stability_archive_path,"r:gz") as archive:
        archive.extractall("third_party")