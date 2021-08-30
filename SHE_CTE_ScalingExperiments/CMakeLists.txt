# Declare this directory as an Elements module
elements_subdir(SHE_CTE_ScalingExperiments)

# Declare Elements module dependencies
elements_depends_on_subdirs(ElementsKernel)
elements_depends_on_subdirs(SHE_CTE)
elements_depends_on_subdirs(SHE_PPT)

# Instruction for Python module installation
elements_install_python_modules()

# Instruction for creating a Python executable
elements_add_python_program(SHE_CTE_SplitFits SHE_CTE_ScalingExperiments.SplitFits)
elements_add_python_program(SHE_CTE_CombineSplitFitsListfile SHE_CTE_ScalingExperiments.CombineSplitFitsListfile)
elements_add_python_program(SHE_CTE_ExtractObjects SHE_CTE_ScalingExperiments.ExtractObjects)
elements_add_python_program(SHE_CTE_MakeBatches SHE_CTE_ScalingExperiments.MakeBatches)
elements_add_python_program(SHE_CTE_ExtractStamps SHE_CTE_ScalingExperiments.ExtractStamps)
elements_add_python_program(SHE_CTE_AnalyseRuntime SHE_CTE_ScalingExperiments.AnalyseRuntime)

# Install the configuration files
elements_install_conf_files()

# Install the auxiliary files
# elements_install_aux_files()
                
