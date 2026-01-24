#pragma once

namespace jzzz
{
    template<typename T>
    void param_class<T>::get_max_cpus()
    {
        const char *slurm_cpus=getenv("SLURM_CPUS_PER_TASK");
        int available_cpus=0;
        if (slurm_cpus!=nullptr)
            available_cpus=atoi(slurm_cpus);
        else
            throw runtime_error("没有获取到每个task可用的cpu数量");
        if (max_cpus==0)
            max_cpus=available_cpus;
        else
            max_cpus=min(max_cpus,available_cpus);
        if ((uint64_t)max_cpus>element_num_to_process)
            max_cpus=element_num_to_process;
        printf("Max cpus= %d\n",max_cpus);
    }
}