# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/smith/CLionProjects/migrate-3.7.2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/smith/CLionProjects/migrate-3.7.2/build

# Include any dependencies generated for this target.
include CMakeFiles/migrate.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/migrate.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/migrate.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/migrate.dir/flags.make

CMakeFiles/migrate.dir/src/main.c.o: CMakeFiles/migrate.dir/flags.make
CMakeFiles/migrate.dir/src/main.c.o: /home/smith/CLionProjects/migrate-3.7.2/src/main.c
CMakeFiles/migrate.dir/src/main.c.o: CMakeFiles/migrate.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/smith/CLionProjects/migrate-3.7.2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/migrate.dir/src/main.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/migrate.dir/src/main.c.o -MF CMakeFiles/migrate.dir/src/main.c.o.d -o CMakeFiles/migrate.dir/src/main.c.o -c /home/smith/CLionProjects/migrate-3.7.2/src/main.c

CMakeFiles/migrate.dir/src/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/migrate.dir/src/main.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/smith/CLionProjects/migrate-3.7.2/src/main.c > CMakeFiles/migrate.dir/src/main.c.i

CMakeFiles/migrate.dir/src/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/migrate.dir/src/main.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/smith/CLionProjects/migrate-3.7.2/src/main.c -o CMakeFiles/migrate.dir/src/main.c.s

# Object files for target migrate
migrate_OBJECTS = \
"CMakeFiles/migrate.dir/src/main.c.o"

# External object files for target migrate
migrate_EXTERNAL_OBJECTS =

migrate: CMakeFiles/migrate.dir/src/main.c.o
migrate: CMakeFiles/migrate.dir/build.make
migrate: src/zlib/libzlib.a
migrate: src/haru/libharu.a
migrate: CMakeFiles/migrate.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/smith/CLionProjects/migrate-3.7.2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable migrate"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/migrate.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/migrate.dir/build: migrate
.PHONY : CMakeFiles/migrate.dir/build

CMakeFiles/migrate.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/migrate.dir/cmake_clean.cmake
.PHONY : CMakeFiles/migrate.dir/clean

CMakeFiles/migrate.dir/depend:
	cd /home/smith/CLionProjects/migrate-3.7.2/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/smith/CLionProjects/migrate-3.7.2 /home/smith/CLionProjects/migrate-3.7.2 /home/smith/CLionProjects/migrate-3.7.2/build /home/smith/CLionProjects/migrate-3.7.2/build /home/smith/CLionProjects/migrate-3.7.2/build/CMakeFiles/migrate.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/migrate.dir/depend
