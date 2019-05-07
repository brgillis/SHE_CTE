# Declare this directory as an Elements module
elements_subdir(SHE_CTE_PipelineUtility)

# Declare Elements module dependencies
elements_depends_on_subdirs(ElementsKernel)
elements_depends_on_subdirs(SHE_CTE)
elements_depends_on_subdirs(SHE_PPT)

# Instruction for Python module installation
elements_install_python_modules()

# Instruction for creating a Python executable
elements_add_python_program(SHE_CTE_ObjectIdSplit SHE_CTE_PipelineUtility.ObjectIdSplit)
elements_add_python_program(SHE_CTE_ShearEstimatesMerge SHE_CTE_PipelineUtility.ShearEstimatesMerge)
elements_add_python_program(SHE_CTE_MakeMockAnalysisData SHE_CTE_PipelineUtility.MakeMockAnalysisData)
elements_add_python_program(SHE_CTE_CleanupBiasMeasurement SHE_CTE_PipelineUtility.CleanupBiasMeasurement)

# Install the configuration files
elements_install_conf_files()

# Install the auxiliary files
# elements_install_aux_files()
                
