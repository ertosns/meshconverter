cmake_minimum_required(VERSION 2.3)

PROJECT(MESHLIB C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CUR  "${CMAKE_CURRENT_SOURCE_DIR}")
set(INCLUDE  "${CUR}/include")
set(SRC  "${CUR}/src")

include_directories(
	"${INCLUDE}"
	"${INCLUDE}/io"
	"${INCLUDE}/converters"
	"${SRC}"
	"${SRC}/IO"
	"${SRC}/converters"
)

file(GLOB_RECURSE HEADER
	"${INCLUDE}/*.h"
	"${INCLUDE}/io/*.h"
	"${INCLUDE}/converters/*.h"
)

file(GLOB_RECURSE SOURCES
	"${SRC}/*.cpp"
	"${SRC}/*.cc"
	"${SRC}/io/*.cpp"
	"${SRC}/converters/*.cpp}")

add_executable(MESHLIB
	${HEADER}
	${SOURCES}
)


