cmake_minimum_required(VERSION 3.10)
project(foma C)

set(foma_VERSION_MAJOR 0)
set(foma_VERSION_MINOR 9)
set(foma_VERSION_PATCH 18)

find_package(FLEX)
find_package(BISON)
find_package(ZLIB REQUIRED)

set(CMAKE_C_STANDARD 99)

if (${BISON_FOUND} AND ${FLEX_FOUND})
	# bison and flex
	bison_target(REGEX "${CMAKE_CURRENT_SOURCE_DIR}/foma/regex.y" "${CMAKE_CURRENT_BINARY_DIR}/regex.c"
		DEFINES_FILE "${CMAKE_CURRENT_BINARY_DIR}/regex.h")
	flex_target(YY "${CMAKE_CURRENT_SOURCE_DIR}/foma/regex.l" "${CMAKE_CURRENT_BINARY_DIR}/lex.yy.c")
	add_flex_bison_dependency(YY REGEX)
	flex_target(LEXC "${CMAKE_CURRENT_SOURCE_DIR}/foma/lexc.l" "${CMAKE_CURRENT_BINARY_DIR}/lex.lexc.c"
		COMPILE_FLAGS "-8 --prefix=lexc")
	flex_target(INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/foma/interface.l" "${CMAKE_CURRENT_BINARY_DIR}/lex.interface.c"
		COMPILE_FLAGS "-8 --prefix=interface")
	flex_target(CMATRIX "${CMAKE_CURRENT_SOURCE_DIR}/foma/cmatrix.l" "${CMAKE_CURRENT_BINARY_DIR}/lex.cmatrix.c"
		COMPILE_FLAGS "-8 --prefix=cmatrix")
endif ()

# source compile
file(GLOB_RECURSE HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/foma" "foma*.h")

add_library(foma OBJECT 
	"${CMAKE_CURRENT_SOURCE_DIR}/foma/apply.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/foma/coaccessible.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/foma/constructions.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/foma/define.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/foma/determinize.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/foma/dynarray.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/foma/extract.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/foma/flags.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/foma/int_stack.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/foma/io.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/foma/lexcread.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/foma/mem.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/foma/minimize.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/foma/reverse.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/foma/sigma.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/foma/spelling.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/foma/stringhash.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/foma/structures.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/foma/topsort.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/foma/trie.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/foma/utf8.c"
	)
target_compile_options(foma PRIVATE -O3 -Wall -D_GNU_SOURCE -fvisibility=hidden -fPIC)
target_include_directories(foma PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/foma)
target_link_libraries(foma PUBLIC ZLIB::ZLIB)

if (${BISON_FOUND} AND ${FLEX_FOUND})
	add_library(foma_script OBJECT
		"${CMAKE_CURRENT_SOURCE_DIR}/foma/rewrite.c"
		${BISON_REGEX_OUTPUT_SOURCE}
		${FLEX_YY_OUTPUTS}
		${FLEX_LEXC_OUTPUTS}
		${FLEX_CMATRIX_OUTPUTS}
		)
	target_compile_options(foma_script PRIVATE -O3 -Wall -D_GNU_SOURCE -fvisibility=hidden -fPIC)
	target_include_directories(foma_script PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/foma)
	target_link_libraries(foma_script PUBLIC ZLIB::ZLIB -lreadline -ltermcap)
endif ()

if (${BISON_FOUND} AND ${FLEX_FOUND})
	add_library(foma_static STATIC $<TARGET_OBJECTS:foma> $<TARGET_OBJECTS:foma_script>)
	target_link_libraries(foma_static PUBLIC ZLIB::ZLIB -lreadline -ltermcap)
else ()
	add_library(foma_static STATIC $<TARGET_OBJECTS:foma>)
	target_link_libraries(foma_static PUBLIC ZLIB::ZLIB)
endif ()
set_target_properties(foma_static PROPERTIES
	OUTPUT_NAME foma
	VERSION 0.9.18
	SOVERSION 0
	PUBLIC_HEADER "${HEADERS}")

if (${BISON_FOUND} AND ${FLEX_FOUND})
	add_library(foma_shared SHARED $<TARGET_OBJECTS:foma> $<TARGET_OBJECTS:foma_script>)
	target_link_libraries(foma_shared PUBLIC ZLIB::ZLIB -lreadline -ltermcap PRIVATE foma)
else ()
	add_library(foma_shared SHARED $<TARGET_OBJECTS:foma>)
	target_link_libraries(foma_shared PUBLIC ZLIB::ZLIB PRIVATE foma)
endif ()
set_target_properties(foma_shared PROPERTIES
	OUTPUT_NAME foma
	VERSION 0.9.18
	SOVERSION 0
	PUBLIC_HEADER "${HEADERS}")

if (${BISON_FOUND} AND ${FLEX_FOUND})
	add_executable(foma_exec
		"${CMAKE_CURRENT_SOURCE_DIR}/foma/foma.c"
		"${CMAKE_CURRENT_SOURCE_DIR}/foma/stack.c"
		"${CMAKE_CURRENT_SOURCE_DIR}/foma/iface.c"
		${FLEX_INTERFACE_OUTPUTS})
	target_include_directories(foma_exec PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/foma)
	target_compile_options(foma_exec PUBLIC -O3 -Wall -D_GNU_SOURCE -fvisibility=hidden)
	target_link_libraries(foma_exec PUBLIC foma ZLIB::ZLIB -lreadline -ltermcap) 
	set_target_properties(foma_exec PROPERTIES OUTPUT_NAME foma)

	add_executable(flookup "${CMAKE_CURRENT_SOURCE_DIR}/foma/flookup.c")
	target_compile_options(flookup PUBLIC -O3 -Wall -D_GNU_SOURCE -fvisibility=hidden)
	target_link_libraries(flookup PUBLIC foma ZLIB::ZLIB -lreadline -ltermcap)

	add_executable(cgflookup "${CMAKE_CURRENT_SOURCE_DIR}/foma/cgflookup.c")
	target_compile_options(cgflookup PUBLIC -O3 -Wall -D_GNU_SOURCE -fvisibility=hidden)
	target_link_libraries(cgflookup PUBLIC foma ZLIB::ZLIB -lreadline -ltermcap)

	install(TARGETS foma flookup cgflookup
		RUNTIME DESTINATION bin)
endif ()

# install
install(TARGETS foma_shared foma_static
	EXPORT foma
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib
	PUBLIC_HEADER DESTINATION include/foma)

include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/foma-config.cmake.in
	${CMAKE_CURRENT_BINARY_DIR}/foma-config.cmake
	INSTALL_DESTINATION lib/cmake/foma)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/foma-config.cmake
	DESTINATION lib/cmake/foma)

install(EXPORT foma 
	NAMESPACE foma::
	DESTINATION lib/cmake/foma)

