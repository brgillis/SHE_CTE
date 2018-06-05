# Declare this directory as an Elements module
elements_subdir(SHE_CTE_PSFFitting)

# Declare Elements module dependencies
elements_depends_on_subdirs(ElementsKernel)
elements_depends_on_subdirs(SHE_CTE)
elements_depends_on_subdirs(SHE_PPT)

# Instruction for Python module installation
elements_install_python_modules()

# Instruction for creating Python executables
elements_add_python_program(SHE_CTE_FitPSFs SHE_CTE_PSFFitting.FitPSFs)
elements_add_python_program(SHE_CTE_ModelPSFs SHE_CTE_PSFFitting.ModelPSFs)

# Install the configuration files
# elements_install_conf_files()

# Install the auxiliary files
elements_install_aux_files()
                
