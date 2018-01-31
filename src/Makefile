SRC  = complex.cxx  iirApply.cxx  iirDesign.cxx  iirFilter.cxx
OBJ  = $(foreach src, ${SRC},  $(subst .cxx,.o,${src}))
LIB = libiirfilter.a

DEBUG = -O3
CXX = g++

all: ${OBJ} ${LIB}


${OBJ}:%.o:%.cxx complex.hh  iirApply.hh  iirDesign.hh  iirFilter.hh
	${CXX} -c -o $@ ${INC} ${DEBUG} -lm $<

${LIB}: ${OBJ}
	ar cr $@ ${OBJ}


clean:
	rm ${OBJ} ${LIB}