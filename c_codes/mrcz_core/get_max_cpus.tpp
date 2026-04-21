#pragma once

namespace mrcz
{
    template<typename T>
    void param_class<T>::get_max_cpus()
    {
        const char *cpus_per_task=getenv("SLURM_CPUS_PER_TASK");
        const char *ntasks=getenv("SLURM_NTASKS");
        int available_cpus=0;
        if(cpus_per_task)
            available_cpus=atoi(cpus_per_task);
        else if(ntasks)
            available_cpus=atoi(ntasks);
        else 
            available_cpus=1;
        if (max_cpus==0)
            max_cpus=available_cpus;
        else
            max_cpus=min(max_cpus,available_cpus);
        if ((uint64_t)max_cpus>element_num_to_process)
            max_cpus=element_num_to_process;
        printf("Max cpus= %d\n",max_cpus);
    }
}