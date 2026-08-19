#pragma once

#include <stdexcept>

#include "mrcz.hpp"

namespace mrcz
{
    template<typename T>
    void validate_common_params(const param_class<T> &params)
    {
        if (params.format!="nekrs" and params.format!="mfem" and params.format!="nektar++")
            throw runtime_error("必须指定数据格式，目前支持: nekrs, mfem, nektar++");
        if (params.data_path.empty())
            throw runtime_error("必须指定数据路径");
        if (params.format=="mfem" and params.mesh_path.empty())
            throw runtime_error("使用mfem格式时必须指定网格文件路径");
        if (params.eb_mode!="ABS" and params.eb_mode!="REL")
            throw runtime_error("误差控制模式必须是ABS或REL");
        if (params.input_eb<=0)
            throw runtime_error("误差值必须为正数");
        if (params.max_cpus<=0)
            throw runtime_error("max_cpus必须是正整数");
        if (params.got_element_range and params.element_begin_index>=params.element_end_index)
            throw runtime_error("element范围必须满足begin<end");
    }
}
