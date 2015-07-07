import os

env = Environment(ENV = os.environ,
                  CXX='clang++',
                  #CXX='g++',
                  CPPPATH=[os.environ['RICH_ROOT']+'/source',
                           os.environ['RICH_ROOT']],
                  LIBPATH=[os.environ['RICH_ROOT'],'.',os.environ['HDF5_LIB_PATH']],
                  LIBS=['rich','hdf5','hdf5_cpp','gfortran'],
                  #LINKFLAGS=' -g -Og ',
                  LINKFLAGS='',
                  #F90FLAGS=' -g -Og -ffpe-trap=invalid,overflow,underflow,zero,denormal -ffpe-summary=all ',
                  F90FLAGS=' -O2 ',
                  #CXXFLAGS='-Weverything -Werror -ferror-limit=1 -Wno-error=padded -g -O0 ')
                  #CXXFLAGS=' -g -Og -Wfatal-errors ')
                  CXXFLAGS=' -Weverything -Werror -ferror-limit=1 -Wno-error=padded -O2 ')
                  
env.Program('rich',Glob('*.cpp'))
