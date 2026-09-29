include_guard(GLOBAL)

function(or_apply_project_options target_name)
    if(NOT TARGET "${target_name}")
        message(FATAL_ERROR "or_apply_project_options: unknown target ${target_name}")
    endif()

    if(OR_ENABLE_WARNINGS)
        if(MSVC)
            target_compile_options("${target_name}" PRIVATE /W4 /permissive-)
        elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
            target_compile_options("${target_name}" PRIVATE
                    -Wall
                    -Wextra
                    -Wpedantic
            )
        endif()
    endif()

    if(OR_ENABLE_ASAN)
        if(MSVC)
            target_compile_options("${target_name}" PRIVATE /fsanitize=address)
        elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
            target_compile_options("${target_name}" PRIVATE
                    -fsanitize=address
                    -fno-omit-frame-pointer
            )
            target_link_options("${target_name}" PRIVATE -fsanitize=address)
        else()
            message(FATAL_ERROR
                    "OR_ENABLE_ASAN is not configured for compiler "
                    "${CMAKE_CXX_COMPILER_ID}")
        endif()
    endif()
endfunction()
