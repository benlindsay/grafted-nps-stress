#CC = g++
CC = mpic++
#CC = ${HOME}/Install/openmpi/bin/mpic++
DFLAGS = -DPAR
CFLAGS = -O3 -w -Wall
INCLUDE = -I${HOME}/Install/fftw3-intel-mpi/include 
LDFLAGS = -limf -lm -lfftw3_mpi -lfftw3 -L${HOME}/Install/fftw3-intel-mpi/lib

SRCS = main.cpp calc_poly_density.cpp homopolymer_discrete.cpp \
       io_utils.cpp integrate_utils.cpp ran2.cpp array_utils.cpp calc_h.cpp \
       fft_wrappers.cpp calc_debye.cpp initialize.cpp \
       read_input.cpp diblock_discrete.cpp 1s_update.cpp \
       Euler_update.cpp brent.cpp simulate.cpp
			 

OBJS = ${SRCS:.cpp=.o}

.cpp.o:
	${CC} ${CFLAGS} ${DFLAGS} -c ${INCLUDE} $<

a.out:	${OBJS}
	${CC} ${CFLAGS} -o $@ ${OBJS} ${LDFLAGS}

clean:
	rm -f *.o a.out
