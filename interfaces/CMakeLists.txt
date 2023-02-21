# Check whether the user has Stimulus on its PATH
find_program(STIMULUS_COMMAND stimulus)

# Add a custom targer that checks functions.yaml and types.yaml
if(STIMULUS_COMMAND)
    add_custom_command(
        OUTPUT test_stimulus_specifications.cpp
        COMMAND ${STIMULUS_COMMAND}
            -l ci:validate
            -f ${CMAKE_CURRENT_SOURCE_DIR}/functions.yaml
            -t ${CMAKE_CURRENT_SOURCE_DIR}/types.yaml
            -o test_stimulus_specifications.cpp
        DEPENDS
            ${CMAKE_CURRENT_SOURCE_DIR}/functions.yaml
            ${CMAKE_CURRENT_SOURCE_DIR}/types.yaml
        COMMENT "Generating C++ checker for Stimulus function and type specifications..."
        USES_TERMINAL
    )
    add_executable(test_stimulus test_stimulus_specifications.cpp)
    target_include_directories(test_stimulus PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/include)

    add_custom_target(
        check_stimulus
        COMMAND test_stimulus
        DEPENDS test_stimulus
        COMMENT "Running C++ checker for Stimulus function and type specifications..."
    )
endif(STIMULUS_COMMAND)
