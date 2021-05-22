FROM registry.gitlab.com/frontistr-commons/frontistr/x86_64-w64-mingw32/openblas:base_hyb AS lib2
RUN curl -L http://www.netlib.org/scalapack/scalapack-2.1.0.tgz | tar zxv && cd scalapack-2.1.0 \
 && sed -e "s/mpif90/x86_64-w64-mingw32-gfortran/g" -e "s/mpicc/x86_64-w64-mingw32-gcc/"  -e "s/ranlib/x86_64-w64-mingw32-ranlib/" -e "s/ar/x86_64-w64-mingw32-ar/" -e "s|-O3|-O3 -I${LIB_ROOT}/include|g" SLmake.inc.example > SLmake.inc \
 && make lib && cp libscalapack.a ${LIB_ROOT}/lib \
 && cd .. && rm -fr scalapack-2.1.0 \
 && curl -L http://mumps.enseeiht.fr/MUMPS_5.4.0.tar.gz | tar zxv && cd MUMPS_5.4.0 \
 && cp Make.inc/Makefile.inc.generic Makefile.inc \
 && make -C src build_mumps_int_def.o build_mumps_int_def \
 && sed \
 -e "s|^CC.*$|CC = ${target}-gcc|"  \
 -e "s|^FC.*$|FC = ${target}-gfortran|"  \
 -e "s|^FL.*$|FL = ${target}-gfortran|" \
 -e "s|^INCPAR.*$|INCPAR = -I${LIB_ROOT}/include|" \
 -e "s|^OPTF.*$|OPTF = -O -fopenmp -DBLR_MT|" \
 -e "s|^OPTC.*$|OPTC = -O -I. -fopenmp|" \
 -e "s|^OPTL.*$|OPTL = -O -fopenmp|" -i Makefile.inc \
 && make RANLIB=${target}-ranlib prerequisites libseqneeded -j \
 && make -C src RANLIB=${target}-ranlib all -j \
 && cp include/*.h ${LIB_ROOT}/include && cp lib/*.a ${LIB_ROOT}/lib \
 && cd .. && rm -fr MUMPS_5.4.0

FROM lib2 AS lib
RUN git clone --depth 1 -b trilinos-release-13-0-1 https://github.com/trilinos/Trilinos.git && cd Trilinos \
 && sed -i -e "s/git.cmd/git/" ./cmake/tribits/core/package_arch/TribitsConstants.cmake \
 && sed -e '1s/^/#include <windows.h>\n/' -e '1s/^/#include <unistd.h>\n/' -i packages/ml/src/Utils/ml_epetra_utils.cpp \
 && mkdir build; cd build \
 && cmake \
  -DCMAKE_TOOLCHAIN_FILE=${LIB_ROOT}/toolchain.cmake -DCMAKE_INSTALL_PREFIX=${LIB_ROOT} \
  -DTPL_ENABLE_MPI=ON \
  -DTrilinos_ENABLE_OpenMP=ON \
  -DCMAKE_CXX_FLAGS_NONE_OVERRIDE=-fopenmp \
  -DBUILD_SHARED_LIBS=OFF -DTPL_ENABLE_DLlib=OFF \
  -DTPL_ENABLE_METIS=ON -DTPL_ENABLE_MUMPS=ON \
  -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF \
  -DTrilinos_ENABLE_TriKota=OFF \
  -DTrilinos_ENABLE_ML=ON \
  -DTrilinos_ENABLE_Zoltan=ON \
  -DTrilinos_ENABLE_Amesos=ON \
  -DBLAS_LIBRARY_NAMES="openblas" -DLAPACK_LIBRARY_NAMES="openblas" \
  -DTrilinos_ENABLE_Fortran=OFF \
  -DHAVE_GCC_ABI_DEMANGLE=1 -DHAVE_TEUCHOS_BLASFLOAT=1 -DHAVE_TEUCHOS_LAPACKLARND=1 \
  -DMPI_C_HEADER_DIR=$LIB_ROOT/include -DMPI_CXX_HEADER_DIR=$LIB_ROOT/include \
  -DMPI_C_ADDITIONAL_INCLUDE_DIRS=$LIB_ROOT/include  -DMPI_CXX_ADDITIONAL_INCLUDE_DIRS=$LIB_ROOT/include \
  -DMPI_C_LIB_NAMES=msmpi -DMPI_CXX_LIB_NAMES=msmpi -DCMAKE_CXX_FLAGS=-fpermissive .. \
 && make -j && make install \
 && cd ../.. && rm -fr Trilinos


FROM lib2 AS lib-trilinos12
RUN git clone --depth 1 -b trilinos-release-12-18-1 https://github.com/trilinos/Trilinos.git \
 && cd Trilinos \
 && sed -i -e "s/git.cmd/git/" ./cmake/tribits/core/package_arch/TribitsConstants.cmake \
 && sed -e '1s/^/#include <windows.h>\n/' -e '1s/^/#include <unistd.h>\n/' -i packages/ml/src/Utils/ml_epetra_utils.cpp \
 && mkdir build; cd build \
 && cmake \
  -DCMAKE_TOOLCHAIN_FILE=${LIB_ROOT}/toolchain.cmake -DCMAKE_INSTALL_PREFIX=${LIB_ROOT} \
  -DTPL_ENABLE_MPI=ON \
  -DTrilinos_ENABLE_OpenMP=ON \
  -DCMAKE_CXX_FLAGS_NONE_OVERRIDE=-fopenmp \
  -DBUILD_SHARED_LIBS=OFF -DTPL_ENABLE_DLlib=OFF \
  -DTPL_ENABLE_METIS=ON -DTPL_ENABLE_MUMPS=ON \
  -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF \
  -DTrilinos_ENABLE_TriKota=OFF \
  -DTrilinos_ENABLE_ML=ON \
  -DTrilinos_ENABLE_Zoltan=ON \
  -DTrilinos_ENABLE_Amesos=ON \
  -DBLAS_LIBRARY_NAMES="openblas" -DLAPACK_LIBRARY_NAMES="openblas" \
  -DTrilinos_ENABLE_Fortran=OFF \
  -DHAVE_GCC_ABI_DEMANGLE=1 -DHAVE_TEUCHOS_BLASFLOAT=1 -DHAVE_TEUCHOS_LAPACKLARND=1 \
  -DMPI_C_HEADER_DIR=$LIB_ROOT/include -DMPI_CXX_HEADER_DIR=$LIB_ROOT/include \
  -DMPI_C_ADDITIONAL_INCLUDE_DIRS=$LIB_ROOT/include  -DMPI_CXX_ADDITIONAL_INCLUDE_DIRS=$LIB_ROOT/include \
  -DMPI_C_LIB_NAMES=msmpi -DMPI_CXX_LIB_NAMES=msmpi -DCMAKE_CXX_FLAGS=-fpermissive .. \
 && make -j && make install \
 && cd ../.. && rm -fr Trilinos

FROM lib2 AS lib-metis4
RUN curl -L http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz | tar zxv && cd metis-4.0.3 \
 && sed -e 's/CC = cc/CC = x86_64-w64-mingw32-gcc/' -i Makefile.in && sed -e 's/COPTIONS = /COPTIONS = -D__VC__/' -i Makefile.in \
 && make -C Lib -j && cp libmetis.a $LIB_ROOT/lib/ \
 && find Lib -name "*.h"|xargs -i cp {} $LIB_ROOT/include/ \
 && cd .. && rm -fr metis-4.0.3 \
 && git clone --depth 1 -b trilinos-release-12-18-1 https://github.com/trilinos/Trilinos.git \
 && cd Trilinos \
 && sed -i -e "s/git.cmd/git/" ./cmake/tribits/core/package_arch/TribitsConstants.cmake \
 && sed -e '1s/^/#include <windows.h>\n/' -e '1s/^/#include <unistd.h>\n/' -i packages/ml/src/Utils/ml_epetra_utils.cpp \
 && mkdir build; cd build \
 && cmake \
  -DCMAKE_TOOLCHAIN_FILE=${LIB_ROOT}/toolchain.cmake -DCMAKE_INSTALL_PREFIX=${LIB_ROOT} \
  -DTPL_ENABLE_MPI=ON \
  -DTrilinos_ENABLE_OpenMP=ON \
  -DCMAKE_CXX_FLAGS_NONE_OVERRIDE=-fopenmp \
  -DBUILD_SHARED_LIBS=OFF -DTPL_ENABLE_DLlib=OFF \
  -DTPL_ENABLE_METIS=ON -DTPL_ENABLE_MUMPS=ON \
  -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF \
  -DTrilinos_ENABLE_TriKota=OFF \
  -DTrilinos_ENABLE_ML=ON \
  -DTrilinos_ENABLE_Zoltan=ON \
  -DTrilinos_ENABLE_Amesos=ON \
  -DBLAS_LIBRARY_NAMES="openblas" -DLAPACK_LIBRARY_NAMES="openblas" \
  -DTrilinos_ENABLE_Fortran=OFF \
  -DHAVE_GCC_ABI_DEMANGLE=1 -DHAVE_TEUCHOS_BLASFLOAT=1 -DHAVE_TEUCHOS_LAPACKLARND=1 \
  -DMPI_C_HEADER_DIR=$LIB_ROOT/include -DMPI_CXX_HEADER_DIR=$LIB_ROOT/include \
  -DMPI_C_ADDITIONAL_INCLUDE_DIRS=$LIB_ROOT/include  -DMPI_CXX_ADDITIONAL_INCLUDE_DIRS=$LIB_ROOT/include \
  -DMPI_C_LIB_NAMES=msmpi -DMPI_CXX_LIB_NAMES=msmpi -DCMAKE_CXX_FLAGS=-fpermissive .. \
 && make -j && make install \
 && cd ../.. && rm -fr Trilinos
