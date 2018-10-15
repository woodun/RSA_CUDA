NVCC = nvcc
NCXXFLAGS = -Xptxas -dlcm=cg --optimize 0
#NVCCFLAGS = -g -G -m64 -arch compute_20

execs = timing_test rsa example

#LIBS += cudart
#LDLIBS = $(addprefix -l, $(LIBS))
#LDFLAGS = -L/usr/local/cuda/lib64/ -lcudart

###
# Files we don't compile with directly, but we want Make to know that we
# depend on them
###
#THER_FILES = cuda_string.h digit.h mpz.h

default: $(execs)

########################################################
#																											 #
# CPU test program																		 #
#																											 #
########################################################
# timing_test: timing_test.cu
# 		$(NVCC) $(NCXXFLAGS) $@.cu -o $@

# example:
# 		$(NVCC) $(NCXXFLAGS) $@.cu -o $@

$(execs): %: %.cu
	$(NVCC) $(NCXXFLAGS) $< -o $@

dump: $(execs)
	for i in $^; do\
		$(NVCC) -c  $$i.cu;\
		cuobjdump -ptx $$i.o > $$i.dump;\
	done

%.c:
%.h:
%.cu:
%.cpp:

clean:
	rm -f $(execs) *.dump *.o
