# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/moroz/project/ipm

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/moroz/project/ipm/build

# Include any dependencies generated for this target.
include CMakeFiles/intersect.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/intersect.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/intersect.dir/flags.make

CMakeFiles/intersect.dir/main.cpp.o: CMakeFiles/intersect.dir/flags.make
CMakeFiles/intersect.dir/main.cpp.o: ../main.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/moroz/project/ipm/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/intersect.dir/main.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/intersect.dir/main.cpp.o -c /home/moroz/project/ipm/main.cpp

CMakeFiles/intersect.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/intersect.dir/main.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/moroz/project/ipm/main.cpp > CMakeFiles/intersect.dir/main.cpp.i

CMakeFiles/intersect.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/intersect.dir/main.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/moroz/project/ipm/main.cpp -o CMakeFiles/intersect.dir/main.cpp.s

CMakeFiles/intersect.dir/main.cpp.o.requires:
.PHONY : CMakeFiles/intersect.dir/main.cpp.o.requires

CMakeFiles/intersect.dir/main.cpp.o.provides: CMakeFiles/intersect.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/intersect.dir/build.make CMakeFiles/intersect.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/intersect.dir/main.cpp.o.provides

CMakeFiles/intersect.dir/main.cpp.o.provides.build: CMakeFiles/intersect.dir/main.cpp.o
.PHONY : CMakeFiles/intersect.dir/main.cpp.o.provides.build

CMakeFiles/intersect.dir/funcs.cpp.o: CMakeFiles/intersect.dir/flags.make
CMakeFiles/intersect.dir/funcs.cpp.o: ../funcs.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/moroz/project/ipm/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/intersect.dir/funcs.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/intersect.dir/funcs.cpp.o -c /home/moroz/project/ipm/funcs.cpp

CMakeFiles/intersect.dir/funcs.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/intersect.dir/funcs.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/moroz/project/ipm/funcs.cpp > CMakeFiles/intersect.dir/funcs.cpp.i

CMakeFiles/intersect.dir/funcs.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/intersect.dir/funcs.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/moroz/project/ipm/funcs.cpp -o CMakeFiles/intersect.dir/funcs.cpp.s

CMakeFiles/intersect.dir/funcs.cpp.o.requires:
.PHONY : CMakeFiles/intersect.dir/funcs.cpp.o.requires

CMakeFiles/intersect.dir/funcs.cpp.o.provides: CMakeFiles/intersect.dir/funcs.cpp.o.requires
	$(MAKE) -f CMakeFiles/intersect.dir/build.make CMakeFiles/intersect.dir/funcs.cpp.o.provides.build
.PHONY : CMakeFiles/intersect.dir/funcs.cpp.o.provides

CMakeFiles/intersect.dir/funcs.cpp.o.provides.build: CMakeFiles/intersect.dir/funcs.cpp.o
.PHONY : CMakeFiles/intersect.dir/funcs.cpp.o.provides.build

CMakeFiles/intersect.dir/vertices.cpp.o: CMakeFiles/intersect.dir/flags.make
CMakeFiles/intersect.dir/vertices.cpp.o: ../vertices.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/moroz/project/ipm/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/intersect.dir/vertices.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/intersect.dir/vertices.cpp.o -c /home/moroz/project/ipm/vertices.cpp

CMakeFiles/intersect.dir/vertices.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/intersect.dir/vertices.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/moroz/project/ipm/vertices.cpp > CMakeFiles/intersect.dir/vertices.cpp.i

CMakeFiles/intersect.dir/vertices.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/intersect.dir/vertices.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/moroz/project/ipm/vertices.cpp -o CMakeFiles/intersect.dir/vertices.cpp.s

CMakeFiles/intersect.dir/vertices.cpp.o.requires:
.PHONY : CMakeFiles/intersect.dir/vertices.cpp.o.requires

CMakeFiles/intersect.dir/vertices.cpp.o.provides: CMakeFiles/intersect.dir/vertices.cpp.o.requires
	$(MAKE) -f CMakeFiles/intersect.dir/build.make CMakeFiles/intersect.dir/vertices.cpp.o.provides.build
.PHONY : CMakeFiles/intersect.dir/vertices.cpp.o.provides

CMakeFiles/intersect.dir/vertices.cpp.o.provides.build: CMakeFiles/intersect.dir/vertices.cpp.o
.PHONY : CMakeFiles/intersect.dir/vertices.cpp.o.provides.build

CMakeFiles/intersect.dir/basepoints.cpp.o: CMakeFiles/intersect.dir/flags.make
CMakeFiles/intersect.dir/basepoints.cpp.o: ../basepoints.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/moroz/project/ipm/build/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/intersect.dir/basepoints.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/intersect.dir/basepoints.cpp.o -c /home/moroz/project/ipm/basepoints.cpp

CMakeFiles/intersect.dir/basepoints.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/intersect.dir/basepoints.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/moroz/project/ipm/basepoints.cpp > CMakeFiles/intersect.dir/basepoints.cpp.i

CMakeFiles/intersect.dir/basepoints.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/intersect.dir/basepoints.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/moroz/project/ipm/basepoints.cpp -o CMakeFiles/intersect.dir/basepoints.cpp.s

CMakeFiles/intersect.dir/basepoints.cpp.o.requires:
.PHONY : CMakeFiles/intersect.dir/basepoints.cpp.o.requires

CMakeFiles/intersect.dir/basepoints.cpp.o.provides: CMakeFiles/intersect.dir/basepoints.cpp.o.requires
	$(MAKE) -f CMakeFiles/intersect.dir/build.make CMakeFiles/intersect.dir/basepoints.cpp.o.provides.build
.PHONY : CMakeFiles/intersect.dir/basepoints.cpp.o.provides

CMakeFiles/intersect.dir/basepoints.cpp.o.provides.build: CMakeFiles/intersect.dir/basepoints.cpp.o
.PHONY : CMakeFiles/intersect.dir/basepoints.cpp.o.provides.build

CMakeFiles/intersect.dir/integral.cpp.o: CMakeFiles/intersect.dir/flags.make
CMakeFiles/intersect.dir/integral.cpp.o: ../integral.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/moroz/project/ipm/build/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/intersect.dir/integral.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/intersect.dir/integral.cpp.o -c /home/moroz/project/ipm/integral.cpp

CMakeFiles/intersect.dir/integral.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/intersect.dir/integral.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/moroz/project/ipm/integral.cpp > CMakeFiles/intersect.dir/integral.cpp.i

CMakeFiles/intersect.dir/integral.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/intersect.dir/integral.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/moroz/project/ipm/integral.cpp -o CMakeFiles/intersect.dir/integral.cpp.s

CMakeFiles/intersect.dir/integral.cpp.o.requires:
.PHONY : CMakeFiles/intersect.dir/integral.cpp.o.requires

CMakeFiles/intersect.dir/integral.cpp.o.provides: CMakeFiles/intersect.dir/integral.cpp.o.requires
	$(MAKE) -f CMakeFiles/intersect.dir/build.make CMakeFiles/intersect.dir/integral.cpp.o.provides.build
.PHONY : CMakeFiles/intersect.dir/integral.cpp.o.provides

CMakeFiles/intersect.dir/integral.cpp.o.provides.build: CMakeFiles/intersect.dir/integral.cpp.o
.PHONY : CMakeFiles/intersect.dir/integral.cpp.o.provides.build

CMakeFiles/intersect.dir/tgen.cpp.o: CMakeFiles/intersect.dir/flags.make
CMakeFiles/intersect.dir/tgen.cpp.o: ../tgen.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/moroz/project/ipm/build/CMakeFiles $(CMAKE_PROGRESS_6)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/intersect.dir/tgen.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/intersect.dir/tgen.cpp.o -c /home/moroz/project/ipm/tgen.cpp

CMakeFiles/intersect.dir/tgen.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/intersect.dir/tgen.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/moroz/project/ipm/tgen.cpp > CMakeFiles/intersect.dir/tgen.cpp.i

CMakeFiles/intersect.dir/tgen.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/intersect.dir/tgen.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/moroz/project/ipm/tgen.cpp -o CMakeFiles/intersect.dir/tgen.cpp.s

CMakeFiles/intersect.dir/tgen.cpp.o.requires:
.PHONY : CMakeFiles/intersect.dir/tgen.cpp.o.requires

CMakeFiles/intersect.dir/tgen.cpp.o.provides: CMakeFiles/intersect.dir/tgen.cpp.o.requires
	$(MAKE) -f CMakeFiles/intersect.dir/build.make CMakeFiles/intersect.dir/tgen.cpp.o.provides.build
.PHONY : CMakeFiles/intersect.dir/tgen.cpp.o.provides

CMakeFiles/intersect.dir/tgen.cpp.o.provides.build: CMakeFiles/intersect.dir/tgen.cpp.o
.PHONY : CMakeFiles/intersect.dir/tgen.cpp.o.provides.build

CMakeFiles/intersect.dir/graph.cpp.o: CMakeFiles/intersect.dir/flags.make
CMakeFiles/intersect.dir/graph.cpp.o: ../graph.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/moroz/project/ipm/build/CMakeFiles $(CMAKE_PROGRESS_7)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/intersect.dir/graph.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/intersect.dir/graph.cpp.o -c /home/moroz/project/ipm/graph.cpp

CMakeFiles/intersect.dir/graph.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/intersect.dir/graph.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/moroz/project/ipm/graph.cpp > CMakeFiles/intersect.dir/graph.cpp.i

CMakeFiles/intersect.dir/graph.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/intersect.dir/graph.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/moroz/project/ipm/graph.cpp -o CMakeFiles/intersect.dir/graph.cpp.s

CMakeFiles/intersect.dir/graph.cpp.o.requires:
.PHONY : CMakeFiles/intersect.dir/graph.cpp.o.requires

CMakeFiles/intersect.dir/graph.cpp.o.provides: CMakeFiles/intersect.dir/graph.cpp.o.requires
	$(MAKE) -f CMakeFiles/intersect.dir/build.make CMakeFiles/intersect.dir/graph.cpp.o.provides.build
.PHONY : CMakeFiles/intersect.dir/graph.cpp.o.provides

CMakeFiles/intersect.dir/graph.cpp.o.provides.build: CMakeFiles/intersect.dir/graph.cpp.o
.PHONY : CMakeFiles/intersect.dir/graph.cpp.o.provides.build

# Object files for target intersect
intersect_OBJECTS = \
"CMakeFiles/intersect.dir/main.cpp.o" \
"CMakeFiles/intersect.dir/funcs.cpp.o" \
"CMakeFiles/intersect.dir/vertices.cpp.o" \
"CMakeFiles/intersect.dir/basepoints.cpp.o" \
"CMakeFiles/intersect.dir/integral.cpp.o" \
"CMakeFiles/intersect.dir/tgen.cpp.o" \
"CMakeFiles/intersect.dir/graph.cpp.o"

# External object files for target intersect
intersect_EXTERNAL_OBJECTS =

intersect: CMakeFiles/intersect.dir/main.cpp.o
intersect: CMakeFiles/intersect.dir/funcs.cpp.o
intersect: CMakeFiles/intersect.dir/vertices.cpp.o
intersect: CMakeFiles/intersect.dir/basepoints.cpp.o
intersect: CMakeFiles/intersect.dir/integral.cpp.o
intersect: CMakeFiles/intersect.dir/tgen.cpp.o
intersect: CMakeFiles/intersect.dir/graph.cpp.o
intersect: CMakeFiles/intersect.dir/build.make
intersect: CMakeFiles/intersect.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable intersect"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/intersect.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/intersect.dir/build: intersect
.PHONY : CMakeFiles/intersect.dir/build

CMakeFiles/intersect.dir/requires: CMakeFiles/intersect.dir/main.cpp.o.requires
CMakeFiles/intersect.dir/requires: CMakeFiles/intersect.dir/funcs.cpp.o.requires
CMakeFiles/intersect.dir/requires: CMakeFiles/intersect.dir/vertices.cpp.o.requires
CMakeFiles/intersect.dir/requires: CMakeFiles/intersect.dir/basepoints.cpp.o.requires
CMakeFiles/intersect.dir/requires: CMakeFiles/intersect.dir/integral.cpp.o.requires
CMakeFiles/intersect.dir/requires: CMakeFiles/intersect.dir/tgen.cpp.o.requires
CMakeFiles/intersect.dir/requires: CMakeFiles/intersect.dir/graph.cpp.o.requires
.PHONY : CMakeFiles/intersect.dir/requires

CMakeFiles/intersect.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/intersect.dir/cmake_clean.cmake
.PHONY : CMakeFiles/intersect.dir/clean

CMakeFiles/intersect.dir/depend:
	cd /home/moroz/project/ipm/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/moroz/project/ipm /home/moroz/project/ipm /home/moroz/project/ipm/build /home/moroz/project/ipm/build /home/moroz/project/ipm/build/CMakeFiles/intersect.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/intersect.dir/depend

