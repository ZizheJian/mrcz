import os
import re

import numpy as np
from paraview import servermanager
from paraview.simple import Slice,TrivialProducer
from vtkmodules.util.numpy_support import numpy_to_vtk,numpy_to_vtkIdTypeArray,vtk_to_numpy
from vtkmodules.util.vtkConstants import VTK_HEXAHEDRON
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import vtkCellArray,vtkUnstructuredGrid
from vtkmodules.vtkFiltersCore import vtkTriangleFilter

from read_nekrs_header import read_nekrs_header

input_path=os.environ.get("PB0_FIELD_FILE","/anvil/projects/x-cis240192/x-zjian1-in-cis240192/element-data/hyuan/pb0.f00001")
field_name=os.environ.get("PB0_FIELD","P").upper()
if field_name not in ("P","T"):raise ValueError("PB0_FIELD目前只支持P或T")
output_path=os.environ.get("PB0_X0_SLICE_DATA",f"results/visualization/pb0_{field_name}_x0_slice.npz")
field_values_path=os.environ.get("PB0_FIELD_VALUES_FILE")
selected_input_path=os.environ.get("PB0_SELECTED_ELEMENTS_INPUT")
selected_output_path=os.environ.get("PB0_SELECTED_ELEMENTS_OUTPUT")
plane_x=float(os.environ.get("PB0_PLANE_X","0"))
selection_margin=float(os.environ.get("PB0_X_MARGIN","1e-6"))
scan_chunk_size=int(os.environ.get("PB0_SCAN_CHUNK","20000"))

def open_coordinates_and_field(path):
    with open(path,"rb") as input_file:
        header=read_nekrs_header(input_file)
        input_file.seek(header.file_element_num*4,os.SEEK_CUR)
        data_offset=input_file.tell()
    dtype=np.dtype(header.np_dtype).newbyteorder(header.endian_mode)
    element_num=header.file_element_num
    element_size=header.element_size
    itemsize=dtype.itemsize
    coordinates=None
    field_values=None
    variable_index=0
    while variable_index<len(header.variables):
        variable=header.variables[variable_index]
        if variable=="X":
            coordinates=np.memmap(path,dtype=dtype,mode="r",offset=data_offset,shape=(element_num,3,element_size))
            data_offset+=element_num*3*element_size*itemsize
        elif variable=="U":
            data_offset+=element_num*3*element_size*itemsize
        elif variable=="P":
            if field_name=="P":field_values=np.memmap(path,dtype=dtype,mode="r",offset=data_offset,shape=(element_num,element_size))
            data_offset+=element_num*element_size*itemsize
        elif variable=="T":
            if field_name=="T":field_values=np.memmap(path,dtype=dtype,mode="r",offset=data_offset,shape=(element_num,element_size))
            data_offset+=element_num*element_size*itemsize
        elif variable=="S":
            scalar_match=re.match(r"S(\d{2})",header.variables[variable_index:])
            if scalar_match is None:raise ValueError(f"无法解析被动标量说明：{header.variables}")
            scalar_num=int(scalar_match.group(1))
            data_offset+=element_num*scalar_num*element_size*itemsize
            variable_index+=2
        else:
            raise ValueError(f"不支持的NekRS变量说明：{variable}")
        variable_index+=1
    if coordinates is None or field_values is None:raise ValueError(f"输入文件必须同时包含坐标X和field {field_name}")
    return header,coordinates,field_values

def select_intersecting_elements(coordinates):
    selected=[]
    for start in range(0,len(coordinates),scan_chunk_size):
        stop=min(start+scan_chunk_size,len(coordinates))
        x=np.asarray(coordinates[start:stop,0,:])
        mask=(np.min(x,axis=1)<=plane_x+selection_margin)&(np.max(x,axis=1)>=plane_x-selection_margin)
        selected.append(np.flatnonzero(mask).astype(np.int64)+start)
        print(f"Scanned elements: {stop}/{len(coordinates)}, selected: {sum(len(ids) for ids in selected)}",flush=True)
    return np.concatenate(selected)

def make_cell_template(shape):
    nx,ny,nz=shape
    def point_id(x,y,z):return x+nx*(y+ny*z)
    cells=[]
    for z in range(nz-1):
        for y in range(ny-1):
            for x in range(nx-1):
                cells.append((point_id(x,y,z),point_id(x+1,y,z),point_id(x+1,y+1,z),point_id(x,y+1,z),point_id(x,y,z+1),point_id(x+1,y,z+1),point_id(x+1,y+1,z+1),point_id(x,y+1,z+1)))
    return np.asarray(cells,dtype=np.int64)

def make_grid(header,coordinates,field_values,selected):
    points_per_element=header.element_size
    points=np.ascontiguousarray(np.asarray(coordinates[selected,:,:]).transpose(0,2,1).reshape(-1,3))
    selected_values=np.ascontiguousarray(np.asarray(field_values[selected,:]).reshape(-1))
    template=make_cell_template(header.file_element_shape)
    cells_per_element=len(template)
    cell_num=len(selected)*cells_per_element
    connectivity=np.empty(cell_num*8,dtype=np.int64)
    element_chunk_size=1024
    for start in range(0,len(selected),element_chunk_size):
        stop=min(start+element_chunk_size,len(selected))
        bases=np.arange(start,stop,dtype=np.int64)*points_per_element
        connectivity[start*cells_per_element*8:stop*cells_per_element*8]=(bases[:,None,None]+template[None,:,:]).reshape(-1)
    offsets=np.arange(0,len(connectivity)+1,8,dtype=np.int64)
    cell_types=np.full(cell_num,VTK_HEXAHEDRON,dtype=np.uint8)

    vtk_points=vtkPoints()
    vtk_points.SetData(numpy_to_vtk(points,deep=0))
    vtk_cells=vtkCellArray()
    vtk_cells.SetData(numpy_to_vtkIdTypeArray(offsets,deep=0),numpy_to_vtkIdTypeArray(connectivity,deep=0))
    vtk_cell_types=numpy_to_vtk(cell_types,deep=0)
    vtk_field=numpy_to_vtk(selected_values,deep=0)
    vtk_field.SetName(field_name)
    grid=vtkUnstructuredGrid()
    grid.SetPoints(vtk_points)
    grid.SetCells(vtk_cell_types,vtk_cells)
    grid.GetPointData().AddArray(vtk_field)
    grid.GetPointData().SetActiveScalars(field_name)
    print(f"Subset grid: elements={len(selected)}, points={len(points)}, cells={cell_num}",flush=True)
    return grid,(points,selected_values,connectivity,offsets,cell_types)

header,coordinates,field_values=open_coordinates_and_field(input_path)
if field_values_path:
    dtype=np.dtype(header.np_dtype).newbyteorder(header.endian_mode)
    expected_size=header.file_element_num*header.element_size*dtype.itemsize
    actual_size=os.path.getsize(field_values_path)
    if actual_size!=expected_size:raise RuntimeError(f"重建字段文件大小无效: expected={expected_size}, actual={actual_size}")
    field_values=np.memmap(field_values_path,dtype=dtype,mode="r",shape=(header.file_element_num,header.element_size))
print(f"Input: field={field_name}, elements={header.file_element_num}, shape={header.file_element_shape}, dtype={header.np_dtype}, values={field_values_path or input_path}",flush=True)
if selected_input_path:
    selected=np.load(selected_input_path)
    if selected.ndim!=1 or selected.dtype.kind not in "iu":raise RuntimeError("selected element索引无效")
    if len(selected)==0 or selected.min()<0 or selected.max()>=header.file_element_num:raise RuntimeError("selected element索引超出范围")
    selected=selected.astype(np.int64,copy=False)
    print(f"Loaded selected elements: {len(selected)}",flush=True)
else:
    selected=select_intersecting_elements(coordinates)
    if selected_output_path:
        os.makedirs(os.path.dirname(selected_output_path),exist_ok=True)
        np.save(selected_output_path,selected)
        print(f"Saved selected elements: {selected_output_path}",flush=True)
if len(selected)==0:raise RuntimeError("没有element与目标平面相交")
grid,backing_arrays=make_grid(header,coordinates,field_values,selected)

producer=TrivialProducer(registrationName="PB0 x-subset")
producer.GetClientSideObject().SetOutput(grid)
producer.UpdatePipeline()
slice_filter=Slice(registrationName="X=0",Input=producer)
slice_filter.SliceType="Plane"
slice_filter.SliceType.Origin=[plane_x,0.0,0.0]
slice_filter.SliceType.Normal=[1.0,0.0,0.0]
slice_filter.UpdatePipeline()

slice_data=servermanager.Fetch(slice_filter)
triangle_filter=vtkTriangleFilter()
triangle_filter.SetInputData(slice_data)
triangle_filter.Update()
triangle_data=triangle_filter.GetOutput()
points=vtk_to_numpy(triangle_data.GetPoints().GetData())
output_values=vtk_to_numpy(triangle_data.GetPointData().GetArray(field_name))
polygons=triangle_data.GetPolys()
offsets=vtk_to_numpy(polygons.GetOffsetsArray())
connectivity=vtk_to_numpy(polygons.GetConnectivityArray())
slice_bounds=slice_filter.GetDataInformation().GetBounds()
slice_value_range=slice_filter.PointData[field_name].GetRange()
coordinate_tolerance=max(selection_margin,1e-10)
if np.max(np.abs(points[:,0]-plane_x))>coordinate_tolerance:raise RuntimeError("导出的切片点不在目标x平面上")
if not np.allclose([output_values.min(),output_values.max()],slice_value_range,rtol=1e-6,atol=1e-12):raise RuntimeError("导出前后的数值范围不一致")
if offsets[0]!=0 or offsets[-1]!=len(connectivity) or np.any(np.diff(offsets)!=3):raise RuntimeError("三角形连接关系无效")
os.makedirs(os.path.dirname(output_path),exist_ok=True)
np.savez_compressed(output_path,points=points,values=output_values,offsets=offsets,connectivity=connectivity)
print(f"Slice output: points={len(points)}, polygons={len(offsets)-1}",flush=True)
print(f"Slice bounds: {slice_bounds}",flush=True)
print(f"{field_name} range on slice: {slice_value_range}",flush=True)
print(f"Saved: {output_path}",flush=True)
