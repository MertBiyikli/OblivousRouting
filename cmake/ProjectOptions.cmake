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

# ============================================================================
# ASan / OR-Tools ABI compatibility
# ============================================================================
#
# Abseil's absl::flat_hash_map (used internally by OR-Tools' MPSolver, e.g.
# MPSolver::variable_name_to_index_) changes its in-memory layout when the
# translation unit that instantiates it is compiled with AddressSanitizer: it
# adds hidden "generation" tracking fields to catch iterator-invalidation bugs
# (see ABSL_SWISSTABLE_ENABLE_GENERATIONS in
# absl/container/internal/raw_hash_set.h, gated on __SANITIZE_ADDRESS__).
#
# Our prebuilt OR-Tools dependency is NOT built with ASan. Any of our own code
# that IS compiled with ASan and reaches into MPSolver/MPObjective internals
# (e.g. via the inline MPSolver::MutableObjective() accessor) therefore
# computes member offsets that no longer match the actual object layout
# produced by the non-ASan library, causing memory corruption/segfaults deep
# inside libortools (observed as SEGVs inside Abseil's
# raw_hash_set::find_or_prepare_insert, called from
# MPObjective::SetCoefficient).
#
# The fix is to compile the small set of translation units that directly
# touch MPSolver/MPObjective/MPVariable internals WITHOUT sanitizer
# instrumentation, so their view of the OR-Tools types' layout matches the
# prebuilt library. This does not reduce ASan coverage anywhere else in the
# codebase; callers of these translation units are unaffected because they
# only ever interact with OR-Tools-backed solvers through opaque
# pointers/virtual calls, never by inspecting MPSolver's layout directly.
function(or_exempt_sources_from_asan)
    if(NOT OR_ENABLE_ASAN)
        return()
    endif()

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        set_source_files_properties(${ARGN} PROPERTIES
                COMPILE_OPTIONS "-fno-sanitize=address"
        )
    endif()
endfunction()
