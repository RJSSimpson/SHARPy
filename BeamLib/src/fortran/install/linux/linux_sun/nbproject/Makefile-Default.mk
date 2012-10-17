#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
CCADMIN=CCadmin
RANLIB=ranlib
CC=cc
CCC=CC
CXX=CC
FC=sunf77

# Macros
PLATFORM=SunStudio_mrbeam-Linux-x86

# Include project Makefile
include linux_sun2-Makefile.mk

# Object Directory
OBJECTDIR=build/Default/${PLATFORM}

# Object Files
OBJECTFILES=

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	cd .. && make -f Makefile

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf:
	cd .. && make -f Makefile clean

# Subprojects
.clean-subprojects:
