################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../lib/VLM.cpp \
../lib/datatypesx.cpp \
../lib/triads.cpp \
../lib/vorticity.cpp \
../lib/wrapper.cpp 

OBJS += \
./lib/VLM.o \
./lib/datatypesx.o \
./lib/triads.o \
./lib/vorticity.o \
./lib/wrapper.o 

CPP_DEPS += \
./lib/VLM.d \
./lib/datatypesx.d \
./lib/triads.d \
./lib/vorticity.d \
./lib/wrapper.d 


# Each subdirectory must supply rules for building sources it contributes
lib/%.o: ../lib/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/home/rjs10/Downloads/Eigen -I/home/rjs10/SharPyProject/UVLMLib/include -O0 -g3 -Wall -c -fmessage-length=0 -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


