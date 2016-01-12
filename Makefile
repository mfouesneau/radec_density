SRCDIR = ./src
OBJDIR = ./build
CXXFLAGS = -I${SRCDIR} -std=c++11 -O3 -Wall -Wextra -pedantic -DNDEBUG
CPP_FILES := $(wildcard $(SRCDIR)/*.cc)
OBJ_FILES := $(addprefix ${OBJDIR}/,$(notdir $(CPP_FILES:.cc=.o)))
TARGETS = radec_density

all: ${TARGETS}


build/%.o: ${SRCDIR}/%.cc
	@mkdir -p $(@D)
	g++ $(CXXFLAGS) -c -o $@ $<

radec_density: $(OBJ_FILES)
	g++ $(CXXFLAGS) -o $@ $^

clean:
	rm -f ${TARGETS}
	rm -rf ${OBJDIR}


