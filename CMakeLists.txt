CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

#---------------------------------------------------------------
# Load macros and functions for Elements-based projects
find_package(ElementsProject)
#---------------------------------------------------------------

# Declare project name and version
elements_project(SHE_CTE 0.3 USE Elements 5.2.2 SHE_PPT 0.2 SHE_MomentsML 1.0 SHE_BFD 0.1)

