cmake_minimum_required(VERSION 3.26)

aux_source_directory(./classes MAIN_CLASSES)
aux_source_directory(./Equations MAIN_EQUATIONS)
aux_source_directory(./Model MAIN_MODEL)
aux_source_directory(./ModelData MAIN_MODEL_DATA)

include_directories(
        ${CMAKE_CURRENT_BINARY_DIR}
        .
        classes
        Equations
        Model
        ModelData
)


add_library(${PROJECT_NAME}_lib STATIC ${MAIN_CLASSES} ${MAIN_EQUATIONS} ${MAIN_MODEL} ${MAIN_MODEL_DATA})
target_link_libraries(${PROJECT_NAME}_lib ${COVDE_LIB} OpenMP::OpenMP_CXX)

#target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})


add_executable(${PROJECT_NAME} main.cpp )
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_lib)