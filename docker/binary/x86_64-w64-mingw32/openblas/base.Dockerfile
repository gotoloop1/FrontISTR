
FROM registry.gitlab.com/frontistr-commons/frontistr/x86_64-w64-mingw32/base AS ser
RUN git clone --depth 1 -b v0.3.15 https://github.com/xianyi/OpenBLAS.git && cd OpenBLAS \
 && CC=${target}-gcc FC=${target}-gfortran RANLIB=${target}-ranlib HOSTCC=gcc make USE_OPENMP=0 BINARY=64 DYNAMIC_ARCH=1 NO_SHARED=1 -j \
 && make PREFIX=${LIB_ROOT} install \
 && cd .. && rm -fr OpenBLAS \
 && curl -L http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz | tar zxv && cd metis-5.1.0 \
 && sed -i -e "/#include <sys\/resource.h>/d" ./GKlib/gk_arch.h && sed -i -e "/extern int gk_getopt/d" -e "/longopts/d" ./GKlib/gk_getopt.h \
 && mkdir buildserial && cd buildserial \
 && cmake  -DCMAKE_TOOLCHAIN_FILE=$LIB_ROOT/toolchain.cmake -DOPENMP=OFF -DGKRAND=ON -DCMAKE_BUILD_TYPE="Release" -DCMAKE_VERBOSE_MAKEFILE=1  -DGKLIB_PATH=../GKlib -DCMAKE_INSTALL_PREFIX="${LIB_ROOT} " .. \
 && make -j && make install \
 && cd ../.. && rm -fr metis-5.1.0

FROM registry.gitlab.com/frontistr-commons/frontistr/x86_64-w64-mingw32/base AS omp
RUN git clone --depth 1 -b v0.3.15 https://github.com/xianyi/OpenBLAS.git && cd OpenBLAS \
 && CC=${target}-gcc FC=${target}-gfortran RANLIB=${target}-ranlib HOSTCC=gcc LDFLAGS=-fopenmp make USE_OPENMP=1 BINARY=64 DYNAMIC_ARCH=1 NO_SHARED=1 -j \
 && make PREFIX=${LIB_ROOT} install \
 && cd .. && rm -fr OpenBLAS \
 && curl -L http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz | tar zxv && cd metis-5.1.0 \
 && sed -i -e "/#include <sys\/resource.h>/d" ./GKlib/gk_arch.h && sed -i -e "/extern int gk_getopt/d" -e "/longopts/d" ./GKlib/gk_getopt.h \
 && mkdir buildserial && cd buildserial \
 && cmake  -DCMAKE_TOOLCHAIN_FILE=$LIB_ROOT/toolchain.cmake -DOPENMP=ON -DGKRAND=ON -DCMAKE_BUILD_TYPE="Release" -DCMAKE_VERBOSE_MAKEFILE=1  -DGKLIB_PATH=../GKlib -DCMAKE_INSTALL_PREFIX="${LIB_ROOT} " .. \
 && make -j && make install \
 && cd ../.. && rm -fr metis-5.1.0

FROM ser AS mpi
RUN curl -L https://www.frontistr.com/files/ci/msmpi_10.1.2.tar.xz | tar Jxv -C ${LIB_ROOT} \
 && cd ${LIB_ROOT}/lib/ && gendef msmpi.dll && x86_64-w64-mingw32-dlltool -d msmpi.def -l libmsmpi.a -D msmpi.dll && rm msmpi.def && cd -

FROM omp AS hyb
RUN curl -L https://www.frontistr.com/files/ci/msmpi_10.1.2.tar.xz | tar Jxv -C ${LIB_ROOT} \
 && cd ${LIB_ROOT}/lib/ && gendef msmpi.dll && x86_64-w64-mingw32-dlltool -d msmpi.def -l libmsmpi.a -D msmpi.dll && rm msmpi.def && cd -

