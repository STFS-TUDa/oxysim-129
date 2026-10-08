# Defines the macro named check_env_var, which takes an environment variable name as input and checks if the value is true or false. if it doesnt exists it is set to false
macro(get_env_var ENV_VAR_NAME)
    # Retrieves the environment variable and converts its value to uppercase for case-insensitive comparison.    
    string(TOUPPER "$ENV{${ENV_VAR_NAME}}" ENV_VAR_VALUE)

    if(ENV_VAR_VALUE STREQUAL "TRUE")
        set(${ENV_VAR_NAME} ON)
    else()
        set(${ENV_VAR_NAME} OFF)
    endif()
    message(STATUS "${ENV_VAR_NAME} : ${${ENV_VAR_NAME}} ")
endmacro()