#This CMake build file automatically builds Chaste and its third party library dependencies as external projects.
cmake_minimum_required (VERSION 2.8.10)
set(CMAKE_USER_MAKE_RULES_OVERRIDE
   ${CMAKE_CURRENT_SOURCE_DIR}/c_flag_overrides.cmake)
set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX
   ${CMAKE_CURRENT_SOURCE_DIR}/cxx_flag_overrides.cmake)
project(BuildChastePlusThirdPartyLibraries)
include(ExternalProject)

#For GUI configs. Change C, and CXX compiler flags dynamically to static, debug build. 
#The cmake CMAKE_USER_MAKE_RULES_OVERRIDEs above takes care of non GUI builds.
foreach(flag_var
        CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
        CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
   if(${flag_var} MATCHES "/MD")
      string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
   endif(${flag_var} MATCHES "/MD")
endforeach(flag_var)

# A bunch of compiler flags. Since we now propagate the compiler and linker flags to all external projects
# this is one place to set all the flags. 
add_definitions(
-MTd 	# => static debug build
-Z7 	# => embed debugging info in library as opposed to using an external .pdb database
-Yu 	# => use precompiled headers: this prevents a linker warning: "H5FDdirect.obj : warning LNK4221: This object 
		# file does not define any previously undefined public symbols, so it will not be used by any link operation 
		# that consumes this library".
-wd4996 # => disable insecure api warnings
-wd4267 # => disable "possible loss of data due to 'narrowing' conversion, e.g. size_t to int"
-wd4290 # => disable warning "C++ exception specification ignored except to indicate a function is not __declspec(nothrow)"
-wd4005 # => 'identifier' : macro redefinition
-wd4018 # => 'expression' : signed/unsigned mismatch
-wd4244 # => 'argument' : conversion from 'type1' to 'type2', possible loss of data
-wd4101 # => 'identifier' : unreferenced local variable
-wd4554 # => possible operator precedence error warning
-wd4305 # => truncation from type1 to type2
-wd4133 # => 'type' : incompatible types - from 'type1' to 'type2'
-wd4661 # => 'identifier' : no suitable definition provided for explicit template instantiation request
)


if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install" CACHE PATH 
  	"Default install directory. Change as required." FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

set(CHASTE_SOURCE_DIR "${CMAKE_SOURCE_DIR}/chaste" CACHE PATH "Specify the path to the root of Chaste sources.")

set(THIRD_PARTY_SOURCE_DIR "${CMAKE_SOURCE_DIR}/third_party_libs" CACHE PATH 
	"Specify the path to the root of the directory containing CMakeLists.txt that builds the third-party libraries.")


#Build the third-party bilbraries as an external project
externalproject_add(ThirdPartyLibraryBuilder
	SOURCE_DIR ${THIRD_PARTY_SOURCE_DIR}
	CMAKE_GENERATOR ${CMAKE_GENERATOR}
	BINARY_DIR ${THIRD_PARTY_SOURCE_DIR}/build
	#INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/third_party_libs
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}/third_party_libs 
		-DCHASTE_SOURCE_DIR:PATH=${CHASTE_SOURCE_DIR} 
		-DBUILD_BOOST:BOOL=ON
		-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
		-DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}	
  	-DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  	-DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
  	-DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
  	-DCMAKE_C_FLAGS_DEBUG:STRING=${CMAKE_C_FLAGS_DEBUG}
  	-DCMAKE_C_FLAGS_MINSIZEREL:STRING=${CMAKE_C_FLAGS_MINSIZEREL}
  	-DCMAKE_C_FLAGS_RELEASE:STRING=${CMAKE_C_FLAGS_RELEASE}
  	-DCMAKE_C_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_C_FLAGS_RELWITHDEBINFO}
  	-DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
  	-DCMAKE_CXX_FLAGS_DEBUG:STRING=${CMAKE_CXX_FLAGS_DEBUG}
  	-DCMAKE_CXX_FLAGS_MINSIZEREL:STRING=${CMAKE_CXX_FLAGS_MINSIZEREL}
  	-DCMAKE_CXX_FLAGS_RELEASE:STRING=${CMAKE_CXX_FLAGS_RELEASE}
  	-DCMAKE_CXX_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_CXX_FLAGS_RELWITHDEBINFO}

)

#Build Chaste as an external project
externalproject_add(ChasteBuilder
	SOURCE_DIR ${CHASTE_SOURCE_DIR}
	CMAKE_GENERATOR ${CMAKE_GENERATOR}
	BINARY_DIR ${CMAKE_BINARY_DIR}/chaste # Potentially shorten build path. Found out that MSVC does not like too deep build paths, which came up while building
										  # the tests. MSVC 2010 has a limitation of around 260 characters, which was bumped up
										  # to 400, I believe, in MSVC 2012 (and hence why I probably did not hit this previously).
	INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/chaste
	DEPENDS ThirdPartyLibraryBuilder
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}/chaste 
		-DCHASTE_BUILD_ROOT_DIR:PATH=${CMAKE_BINARY_DIR}/chaste
		-DTHIRD_PARTY_LIBS:PATH=${CMAKE_INSTALL_PREFIX}
		-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
		-DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  	-DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  	-DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
  	-DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
  	-DCMAKE_C_FLAGS_DEBUG:STRING=${CMAKE_C_FLAGS_DEBUG}
  	-DCMAKE_C_FLAGS_MINSIZEREL:STRING=${CMAKE_C_FLAGS_MINSIZEREL}
  	-DCMAKE_C_FLAGS_RELEASE:STRING=${CMAKE_C_FLAGS_RELEASE}
  	-DCMAKE_C_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_C_FLAGS_RELWITHDEBINFO}
  	-DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
  	-DCMAKE_CXX_FLAGS_DEBUG:STRING=${CMAKE_CXX_FLAGS_DEBUG}
  	-DCMAKE_CXX_FLAGS_MINSIZEREL:STRING=${CMAKE_CXX_FLAGS_MINSIZEREL}
  	-DCMAKE_CXX_FLAGS_RELEASE:STRING=${CMAKE_CXX_FLAGS_RELEASE}
  	-DCMAKE_CXX_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_CXX_FLAGS_RELWITHDEBINFO}
)

