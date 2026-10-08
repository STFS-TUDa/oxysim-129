# Set the name of your test executable. You can change the name as needed.
set(TEST_NAME testFlameletThermo)

# Add the source file for the test executable.
# Replace ${TEST_NAME}.C with the actual name of your test source file if it's different.
add_executable(${TEST_NAME} 
    ${TEST_NAME}.C
)

# Link the necessary libraries to the test executable, including Catch2.
# Catch2::Catch2WithMain is used to provide the main() function for the tests.
# Add any additional libraries your test requires in this section.
target_link_libraries(${TEST_NAME}  
    PUBLIC 
        Catch2::Catch2WithMain
        # Add any other required libraries here
        flameletThermo
)

# Register the test with CTest.
# The COMMAND can either run the test executable directly or call a script that sets up the test environment.
# If you are using a script, replace `${CMAKE_CURRENT_SOURCE_DIR}/Alltest` with the path to your script.
add_test(NAME ${TEST_NAME} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/Alltest)

# If you want to run the test executable directly instead of a script, uncomment the following line.
# add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})

# Set test properties like the working directory or additional arguments.
# WORKING_DIRECTORY sets where the test will be run. You can modify this if needed.
set_tests_properties(
    ${TEST_NAME}
    PROPERTIES
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)