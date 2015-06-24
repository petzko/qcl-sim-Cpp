SRCDIR=src/
BINDIR=build/
INCS=-I./include -I$(HOME)/libraries/gsl-1.9/install/include
LIBS= -L./lib 


# Include files
SOURCES= matrix.cpp main.cpp
#headers
DEPS= datastructures.hpp 
HEADERS=$(patsubst %.hpp,./include/%.hpp,$(DEPS))


# Compiler
CC=g++
CFLAGS=$(INCS) -Wall -g  

# Linker flags
LDFLAGS= $(LIBS)  -lgsl -lgslcblas 

OBJECTS=$(SOURCES:.cpp=.o)
BINOBJECTS=$(patsubst %.o,$(BINDIR)%.o,$(OBJECTS))
EXECUTABLE=start


LINK=$(CC) $(BINOBJECTS) -o lhs $(LDFLAGS)
COMPILE=$(CC) -c rhs -o $(BINDIR)$@ $(CFLAGS)

all: $(EXECUTABLE)


$(EXECUTABLE): $(OBJECTS)
	@echo  "Link command: "	$(CC) $(BINOBJECTS) -o $@ $(LDFLAGS)
	$(CC) $(BINOBJECTS) -o $@ $(LDFLAGS)

clean:
	rm -f $(BINOBJECTS) $(EXECUTABLE)

$(OBJECTS): %.o : $(SRCDIR)%.cpp $(HEADERS)
	@echo "Headers" = $(HEADERS)
	@echo "Compile command:" $(CC) -c $< -o $(BINDIR)$@ $(CFLAGS)
	$(CC) -c $< $(CFLAGS) -o $(BINDIR)$@
