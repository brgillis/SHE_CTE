CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

#---------------------------------------------------------------
# Load macros and functions for Elements-based projects
find_package(ElementsProject)
#---------------------------------------------------------------

# Declare project name and version
elements_project(SHE_CTE 0.10.6 USE Elements 5.8 SHE_PPT 1.4.1 SHE_BFD 0.8.0 SHE_LensMC 2.10.1 SHE_MomentsML 1.10.1)
