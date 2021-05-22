FROM registry.gitlab.com/frontistr-commons/frontistr/x86_64-w64-mingw32/base AS base
RUN curl -L https://www.frontistr.com/files/ci/mkl_2021.2.0.tar.xz | tar Jxv -C ${LIB_ROOT}

FROM base AS ser
RUN curl -L http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz | tar zxv && cd metis-5.1.0 \
 && sed -i -e "/#include <sys\/resource.h>/d" ./GKlib/gk_arch.h && sed -i -e "/extern int gk_getopt/d" -e "/longopts/d" ./GKlib/gk_getopt.h \
 && mkdir buildserial && cd buildserial \
 && cmake  -DCMAKE_TOOLCHAIN_FILE=$LIB_ROOT/toolchain.cmake -DOPENMP=OFF -DGKRAND=ON -DCMAKE_BUILD_TYPE="Release" -DCMAKE_VERBOSE_MAKEFILE=1  -DGKLIB_PATH=../GKlib -DCMAKE_INSTALL_PREFIX="${LIB_ROOT} " .. \
 && make -j && make install \
 && cd ../.. && rm -fr metis-5.1.0

FROM base AS omp
RUN curl -L http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz | tar zxv && cd metis-5.1.0 \
 && sed -i -e "/#include <sys\/resource.h>/d" ./GKlib/gk_arch.h && sed -i -e "/extern int gk_getopt/d" -e "/longopts/d" ./GKlib/gk_getopt.h \
 && mkdir buildserial && cd buildserial \
 && cmake  -DCMAKE_TOOLCHAIN_FILE=$LIB_ROOT/toolchain.cmake -DOPENMP=ON -DGKRAND=ON -DCMAKE_BUILD_TYPE="Release" -DCMAKE_VERBOSE_MAKEFILE=1  -DGKLIB_PATH=../GKlib -DCMAKE_INSTALL_PREFIX="${LIB_ROOT} " .. \
 && make -j && make install \
 && cd ../.. && rm -fr metis-5.1.0

FROM ser AS mpi
RUN curl -L https://www.frontistr.com/files/ci/impi_2021.2.0.tar.xz | tar Jxv -C ${LIB_ROOT}

FROM omp AS hyb
RUN curl -L https://www.frontistr.com/files/ci/impi_2021.2.0.tar.xz | tar Jxv -C ${LIB_ROOT}