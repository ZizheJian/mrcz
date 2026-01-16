import os

# method可选:
#   "my"
#   "sz3-0": 把数据展成1维的，然后用SZ3的1维压缩
#   "sz3-1": 对每个element执行sz3的3维压缩，在block内共用huffman tree
#   "sz3-2": 对每个element执行sz3的3维压缩，考虑到点的间隔不一致，在block内共用huffman tree
#   "sz2": 对每个element执行sz2的3维压缩，在block内共用huffman tree
#   "zfp": 对每个element执行zfp的3维压缩
#   "sperr": 对每个element执行sperr的3维压缩

method="my"

data_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000"
abs_eb=100

if method=="my":
    if os.path.exists(f"{data_path}.header.yaml"):
        try:
            os.chdir("build")
            os.system("cmake -DCMAKE_INSTALL_PREFIX=../install ../c_codes")
            os.system("make")
            os.system("make install")
            os.chdir("..")
            # os.system(f"./build/jzzz -i {data_path}.header.yaml -M ABS {abs_eb} -blocks [0]")
            os.system(f"./build/jzzz -i {data_path}.header.yaml -M ABS {abs_eb}")
        except Exception as e:
            exit(0)
    else:
        raise NotImplementedError("这里要自动生成header文件和block文件，目前还没写")

