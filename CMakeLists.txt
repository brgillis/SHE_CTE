CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

#---------------------------------------------------------------
# Load macros and functions for Elements-based projects
find_package(ElementsProject)
#---------------------------------------------------------------

# Declare project name, version, and dependencies.
# We only explicitly depend on direct dependencies, and inherit from them dependencies on SHE_PPT and Elements
elements_project(SHE_CTE 8.1 USE SHE_LensMC 3.2) # SHE_BFD 0.9 SHE_MomentsML 1.9)
