# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
  
#Create a variable called linalg_SOURCES containing all .cpp files:

set(miscellaneous_Src 
    src/miscellaneous/allocateSharedVector.c
    src/miscellaneous/allocateSingleSharedVector.c
    PARENT_SCOPE
    )

