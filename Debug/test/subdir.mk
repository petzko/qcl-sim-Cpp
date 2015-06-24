################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../test/bytewisecopytest.cpp \
../test/typeidtest.cpp 

OBJS += \
./test/bytewisecopytest.o \
./test/typeidtest.o 

CPP_DEPS += \
./test/bytewisecopytest.d \
./test/typeidtest.d 


# Each subdirectory must supply rules for building sources it contributes
test/%.o: ../test/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/home/petzko/libraries/gsl-1.9/install/include -I"/home/petzko/workspace/qcl-sim-cpp/include" -I"/home/petzko/workspace/qcl-sim-cpp/include" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


