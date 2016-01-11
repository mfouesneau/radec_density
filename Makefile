SRCDIR = ./src
OBJDIR = ./obj
CXXFLAGS = -I${SRCDIR} -std=c++11 -O3 -Wall -Wextra -pedantic -DNDEBUG
CPP_FILES := $(wildcard $(SRCDIR)/*.cc)
OBJ_FILES := $(addprefix ${OBJDIR}/,$(notdir $(CPP_FILES:.cc=.o)))

all: radec_density

obj/%.o: ${SRCDIR}/%.cc
	   g++ $(CXXFLAGS) -c -o $@ $<

radec_density: $(OBJ_FILES)
	g++ $(CXXFLAGS) -o $@ $^\

clean:
	rm -f radec_density
	rm -f ${OBJDIR}/*.o

