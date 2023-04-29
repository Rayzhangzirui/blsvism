# check plate form, on Mac, compile wiht -Ofast
HOSTNAME = $(shell echo $$HOSTNAME)
#local computer
$(info	HOSTNAME is $(HOSTNAME))

CPPC = g++ -std=c++17 -MMD

ifeq ($(HOSTNAME),MacBook-Air-Ray.local)
omp = 1
endif

ifeq ($(findstring math,$(HOSTNAME)),math)
omp = 0
ocl = 1
OMP = -Xpreprocessor -fopenmp -lomp
PROJECTDIR = /Users/zzirui/c25
endif

ifeq ($(shell hostname | grep tscc | wc -l), 1)
omp = 1
endif

ifeq ($(HOSTNAME),ccom-boom-login.local)
omp = 1
OMP = -fopenmp
endif

# debug flag
ifeq ($(debug), 1)
	CCFLAGS = -g -DDEBUG
else
	CCFLAGS = -Ofast
endif

PROJECTDIR="/Users/Ray/project/c25/vismsrc"
CCFLAGS+=-DPROJECTDIR=\"$(PROJECTDIR)\"

ifeq ($(cpu),1)
	omp=0
	ocl=0
endif

ifeq ($(omp),1)
$(info use omp)
	OPT = -DOPENMP
	SUFFIX = omp
	DEP += ompvism.o
	LIBS += $(OMP)
else ifeq ($(ocl),1)
$(info use ocl)
	OPT = -DOPENCL
	DEP += oclstart.o vismgpu.o
	LIBS += -lm -framework OpenCL
	SUFFIX = ocl
else
$(info use cpu)
	SUFFIX = cpu
endif

DEP += vism_addon.o Solute.o globals.o surf.o cfa.o

DEP += heap.o kernel.o vism.o 

# %.o : %.cpp
# 	$(CPPC) $(CCFLAGS) $(OPT) $(LIBS) -c $< 


.SECONDEXPANSION:
%.o : %.cpp $$(wildcard $$*.h)
	$(CPPC) $(CCFLAGS)  $(OPT) -c $< 

test: test.o $(DEP)
	$(CPPC) $(CCFLAGS) $(OPT) $(LIBS) -lgtest $^ -o $@

test_vism: test_vism.o $(DEP)
	$(CPPC) $(CCFLAGS) $(OPT) $(LIBS) -lgtest -lfmt $^ -o $@

kerneltest: kerneltest.o kernel2.o surf.o
	$(CPPC) $(CCFLAGS) $(OPT) $(LIBS) -lfmt $^ -o $@

atomtest: atomtest.o $(DEP)
	$(CPPC) $(CCFLAGS) $(OPT) $(LIBS) -lfmt $^ -o $@

atom2test: atom2test.o $(DEP)
	$(CPPC) $(CCFLAGS) $(OPT) $(LIBS) -lfmt $^ -o $@

test_kernel: test_kernel.o globals.o kernel.o surf.o
	$(CPPC) $(CCFLAGS) $(OPT) $(LIBS) -lgtest -lfmt  $^ -o $@	 

test_kernel2: test_kernel2.o globals.o kernel2.o surf.o
	$(CPPC) $(CCFLAGS) $(OPT) $(LIBS) -lgtest -lfmt  $^ -o $@	 

test_heap: test_heap.o
	$(CPPC) $(CCFLAGS) $(OPT) $(LIBS) -lgtest -lfmt  $^ -o $@	 

vismtest: vismtest.o $(DEP)
	$(CPPC) $(CCFLAGS) $(OPT) $(LIBS) -lfmt $^ -o $@


test_fft: test_fft.o $(DEP) fftconv.o
	$(CPPC) $(CCFLAGS) $(OPT) $(LIBS)  -lgtest -lfmt -lfftw3 $^ -o $@

vismfft: vismfft.o $(DEP) fftconv.o
	$(CPPC) $(CCFLAGS) $(OPT) $(LIBS) -lfmt -lfftw3 $^ -o $@

all: atomtest atom2test vismtest vismfft

.PHONY : clean
clean :
	rm -rf *.dSYM *.o

cleantest:
	rm $(UTIL) $(TESTS)
