# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
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

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src

# Utility rule file for boost.

# Include the progress variables for this target.
include CMakeFiles/boost.dir/progress.make

CMakeFiles/boost: CMakeFiles/boost-complete


CMakeFiles/boost-complete: boost-prefix/src/boost-stamp/boost-install
CMakeFiles/boost-complete: boost-prefix/src/boost-stamp/boost-mkdir
CMakeFiles/boost-complete: boost-prefix/src/boost-stamp/boost-download
CMakeFiles/boost-complete: boost-prefix/src/boost-stamp/boost-update
CMakeFiles/boost-complete: boost-prefix/src/boost-stamp/boost-patch
CMakeFiles/boost-complete: boost-prefix/src/boost-stamp/boost-configure
CMakeFiles/boost-complete: boost-prefix/src/boost-stamp/boost-build
CMakeFiles/boost-complete: boost-prefix/src/boost-stamp/boost-install
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Completed 'boost'"
	/usr/bin/cmake -E make_directory /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/CMakeFiles
	/usr/bin/cmake -E touch /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/CMakeFiles/boost-complete
	/usr/bin/cmake -E touch /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src/boost-stamp/boost-done

boost-prefix/src/boost-stamp/boost-install: boost-prefix/src/boost-stamp/boost-build
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Performing install step for 'boost'"
	cd /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src/boost && ./b2 variant=release link=static runtime-link=static threading=multi "cxxflags= -std=c++17 -pedantic -O3 " install
	cd /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src/boost && /usr/bin/cmake -E touch /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src/boost-stamp/boost-install

boost-prefix/src/boost-stamp/boost-mkdir:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Creating directories for 'boost'"
	/usr/bin/cmake -E make_directory /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src/boost
	/usr/bin/cmake -E make_directory /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src/boost
	/usr/bin/cmake -E make_directory /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix
	/usr/bin/cmake -E make_directory /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/tmp
	/usr/bin/cmake -E make_directory /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src/boost-stamp
	/usr/bin/cmake -E make_directory /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src
	/usr/bin/cmake -E make_directory /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src/boost-stamp
	/usr/bin/cmake -E touch /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src/boost-stamp/boost-mkdir

boost-prefix/src/boost-stamp/boost-download: boost-prefix/src/boost-stamp/boost-urlinfo.txt
boost-prefix/src/boost-stamp/boost-download: boost-prefix/src/boost-stamp/boost-mkdir
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Performing download step (download, verify and extract) for 'boost'"
	cd /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src && /usr/bin/cmake -P /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src/boost-stamp/download-boost.cmake
	cd /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src && /usr/bin/cmake -P /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src/boost-stamp/verify-boost.cmake
	cd /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src && /usr/bin/cmake -P /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src/boost-stamp/extract-boost.cmake
	cd /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src && /usr/bin/cmake -E touch /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src/boost-stamp/boost-download

boost-prefix/src/boost-stamp/boost-update: boost-prefix/src/boost-stamp/boost-download
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "No update step for 'boost'"
	/usr/bin/cmake -E echo_append
	/usr/bin/cmake -E touch /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src/boost-stamp/boost-update

boost-prefix/src/boost-stamp/boost-patch: boost-prefix/src/boost-stamp/boost-download
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "No patch step for 'boost'"
	/usr/bin/cmake -E echo_append
	/usr/bin/cmake -E touch /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src/boost-stamp/boost-patch

boost-prefix/src/boost-stamp/boost-configure: boost-prefix/tmp/boost-cfgcmd.txt
boost-prefix/src/boost-stamp/boost-configure: boost-prefix/src/boost-stamp/boost-update
boost-prefix/src/boost-stamp/boost-configure: boost-prefix/src/boost-stamp/boost-patch
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Performing configure step for 'boost'"
	cd /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src/boost && CC=/usr/bin/cc CXX=/usr/bin/c++ ./bootstrap.sh --with-toolset=gcc --prefix=/home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost --with-libraries=program_options
	cd /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src/boost && /usr/bin/cmake -E touch /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src/boost-stamp/boost-configure

boost-prefix/src/boost-stamp/boost-build: boost-prefix/src/boost-stamp/boost-configure
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Performing build step for 'boost'"
	cd /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src/boost && ./b2 variant=release link=static runtime-link=static threading=multi "cxxflags= -std=c++17 -pedantic -O3 "
	cd /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src/boost && /usr/bin/cmake -E touch /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/boost-prefix/src/boost-stamp/boost-build

boost: CMakeFiles/boost
boost: CMakeFiles/boost-complete
boost: boost-prefix/src/boost-stamp/boost-install
boost: boost-prefix/src/boost-stamp/boost-mkdir
boost: boost-prefix/src/boost-stamp/boost-download
boost: boost-prefix/src/boost-stamp/boost-update
boost: boost-prefix/src/boost-stamp/boost-patch
boost: boost-prefix/src/boost-stamp/boost-configure
boost: boost-prefix/src/boost-stamp/boost-build
boost: CMakeFiles/boost.dir/build.make

.PHONY : boost

# Rule to build all files generated by this target.
CMakeFiles/boost.dir/build: boost

.PHONY : CMakeFiles/boost.dir/build

CMakeFiles/boost.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/boost.dir/cmake_clean.cmake
.PHONY : CMakeFiles/boost.dir/clean

CMakeFiles/boost.dir/depend:
	cd /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src /home/mehdi/Nextcloud/recherche/code_buss/temp_betweenness_1.0/src/CMakeFiles/boost.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/boost.dir/depend

