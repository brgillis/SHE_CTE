CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

#---------------------------------------------------------------
# Load macros and functions for Elements-based projects
find_package(ElementsProject)
#---------------------------------------------------------------

# Declare project name and version
elements_project(SHE_CTE 0.8.10 USE Elements 5.2.2 SHE_PPT 1.2.2 SHE_BFD 0.6.2 SHE_LensMC 2.8.8 SHE_MomentsML 1.8.5)
