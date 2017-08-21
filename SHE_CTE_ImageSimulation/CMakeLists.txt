# Declare this directory as an Elements module
elements_subdir(SHE_CTE_ImageSimulation)

# Declare Elements module dependencies
elements_depends_on_subdirs(ElementsKernel)
elements_depends_on_subdirs(SHE_GST_GalaxyImageGeneration)
elements_depends_on_subdirs(SHE_GST_IceBRGpy)

# Instruction for Python module installation
elements_install_python_modules()

# Instruction for creating a Python executable
elements_add_python_program(SHE_CTE_SimulateImages SHE_CTE_ImageSimulation.SimulateImages)

# Install the configuration files
# elements_install_conf_files()

# Install the auxiliary files
# elements_install_aux_files()
                
