CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

#---------------------------------------------------------------
# Load macros and functions for Elements-based projects
find_package(ElementsProject)
#---------------------------------------------------------------

# Declare project name and version
elements_project(SHE_CTE 0.9 USE Elements 5.8 SHE_PPT 1.3 SHE_BFD 0.9 SHE_LensMC 2.11 SHE_MomentsML 1.9)
