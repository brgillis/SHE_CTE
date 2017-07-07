# Declare this directory as an Elements module
elements_subdir(SHE_SBM_BiasMeasurement)

# Declare Elements module dependencies
elements_depends_on_subdirs(ElementsKernel)
elements_depends_on_subdirs(icebrgpy)
elements_depends_on_subdirs(SHE_SIM_galaxy_image_generation)
elements_depends_on_subdirs(SHE_SBM_ShearEstimation)

# Instruction for Python module installation
elements_install_python_modules()

# Instruction for creating a Python executable
elements_add_python_program(SHE_SBM_MeasureBias SHE_SBM_BiasMeasurement.MeasureBias)

# Install the configuration files
# elements_install_conf_files()

# Install the auxiliary files
# elements_install_aux_files()
                
