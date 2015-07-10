################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/main.cpp \
../src/matrix.cpp \
../src/matrixview.cpp 

OBJS += \
./src/main.o \
./src/matrix.o \
./src/matrixview.o 

CPP_DEPS += \
./src/main.d \
./src/matrix.d \
./src/matrixview.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/home/petzko/libraries/gsl-1.9/install/include -I"/home/petzko/workspace/qcl-sim-cpp/include" -I"/home/petzko/workspace/qcl-sim-cpp/include" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


