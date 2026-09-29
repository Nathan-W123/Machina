# Write the git revision of SOURCE_DIR into the header OUTPUT as
# SPARLAB_GIT_SHA ("unknown" outside a git checkout; "-dirty" appended when
# tracked files have uncommitted changes). The header is rewritten only when
# the revision changes, so an unchanged tree rebuilds nothing.
set(sha "unknown")
find_package(Git QUIET)
if(GIT_FOUND)
  execute_process(COMMAND ${GIT_EXECUTABLE} -C ${SOURCE_DIR} rev-parse --short=12 HEAD
                  OUTPUT_VARIABLE head OUTPUT_STRIP_TRAILING_WHITESPACE
                  RESULT_VARIABLE status ERROR_QUIET)
  if(status EQUAL 0 AND head)
    set(sha "${head}")
    execute_process(COMMAND ${GIT_EXECUTABLE} -C ${SOURCE_DIR} diff --quiet HEAD --
                    RESULT_VARIABLE dirty ERROR_QUIET)
    if(NOT dirty EQUAL 0)
      set(sha "${sha}-dirty")
    endif()
  endif()
endif()
set(content "#define SPARLAB_GIT_SHA \"${sha}\"\n")
if(EXISTS ${OUTPUT})
  file(READ ${OUTPUT} previous)
endif()
if(NOT "${previous}" STREQUAL "${content}")
  file(WRITE ${OUTPUT} "${content}")
endif()
