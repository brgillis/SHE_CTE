# Declare this directory as an Elements module
elements_subdir(SHE_CTE_BiasMeasurement)

# Declare Elements module dependencies
elements_depends_on_subdirs(ElementsKernel)
elements_depends_on_subdirs(SHE_CTE_PipelineUtility)
elements_depends_on_subdirs(SHE_CTE_ShearEstimation)
elements_depends_on_subdirs(SHE_PPT)
elements_depends_on_subdirs(SHE_BFD_CalculateMoments)

# Instruction for Python module installation
elements_install_python_modules()

# Instruction for creating a Python executable
elements_add_python_program(SHE_CTE_MeasureBias SHE_CTE_BiasMeasurement.MeasureBias)
elements_add_python_program(SHE_CTE_MeasureStatistics SHE_CTE_BiasMeasurement.MeasureStatistics)
elements_add_python_program(SHE_CTE_PlotBias SHE_CTE_BiasMeasurement.PlotBias)
elements_add_python_program(SHE_CTE_PlotPsfSensitivity SHE_CTE_BiasMeasurement.PlotPsfSensitivity)
elements_add_python_program(SHE_CTE_PrintBias SHE_CTE_BiasMeasurement.PrintBias)

# Install the configuration files
# elements_install_conf_files()

# Install the auxiliary files
# elements_install_aux_files()
                
