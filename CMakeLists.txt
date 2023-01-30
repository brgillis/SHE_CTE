CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

#---------------------------------------------------------------
# Load macros and functions for Elements-based projects
find_package(ElementsProject)
#---------------------------------------------------------------

# Declare project name, version, and dependencies.
# We only explicitly depend on direct dependencies, and inherit from them dependencies on SHE_PPT and Elements

elements_project(SHE_CTE 9.2 USE SHE_PPT 9.3.1)
