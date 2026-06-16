CXX = g++



# ROOT
#ROOT_INC  = $(shell root-config --cflags) # liga se for usar ROOT
#ROOT_LIB  = $(shell root-config --libs) # liga se for usar ROOT

#LHAPDF
LHAPDF_CFLAGS = $(shell lhapdf-config --cflags)
LHAPDF_LIBS   = $(shell lhapdf-config --libs)

CXXFLAGS = -O3 -std=c++20 -fopenmp -Wno-deprecated-declarations \
           $(LHAPDF_CFLAGS)
#           $(ROOT_INC)

LDFLAGS = $(PYTHON_LIB) $(LHAPDF_LIBS) $(BOOST_LIBS) #$(ROOT_LIB) 

SRCS = libraries/mantysaari/dipoleamplitude.cpp \
       libraries/mantysaari/dglap_cpp/AlphaStrong.cpp \
       libraries/mantysaari/dglap_cpp/EvolutionLO_nocoupling.cpp \
       main.cpp \
       other/integration.cpp \
       other/utils.cpp \
       other/ctes.cpp \
       dipole_amplitudes/GBW.cpp \
       dipole_amplitudes/IIM.cpp \
       dipole_amplitudes/bCGC.cpp \
       dipole_amplitudes/ipsat.cpp \
       dipole_amplitudes/LHAPDF.cpp \
       other/correcs.cpp \
       calculations/wavefunctions.cpp \
       calculations/nuclear.cpp \
       calculations/sigma.cpp \
       exec/CSV.cpp \

       #dipole_amplitudes/IIM.cpp #
       #dipole_amplitudes/bCGC.cpp #
       #dipole_amplitudes/DGLAP.cpp #

       #ResolveDGLAP.cpp \#
#       minuit.cpp # liga se for fazer ajuste de parâmetros

TARGET = main

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET)
