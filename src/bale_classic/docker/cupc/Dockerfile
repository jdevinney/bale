
FROM ubuntu:18.04

MAINTAINER bale <bale@super.org>

# Build-arg: NPROCS for parallel make
ARG NPROCS=2

ENV PREFIX=/usr/local
RUN mkdir -p "${PREFIX}"
RUN mkdir -p "${PREFIX}/src"


# This will prevent questions from being asked during the install
ENV DEBIAN_FRONTEND noninteractive

# Get certs so HTTPS requests (e.g., git clone) works
# do we need to do this?
RUN apt-get -y update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
    && apt-get clean

# Install various useful/prereq packages 
RUN apt-get -y update \
    && apt-get install -y --no-install-recommends \
        autoconf \
        automake \
	cmake \
        g++ \
        gcc \
        gfortran \
        git \
	libnuma-dev \
        libtool \
        m4 \
	make \
	openssh-client \
        patch \
	pkg-config \
	python \
        wget \
        zlib1g-dev \
    && apt-get clean

#ENV CUPC_URL=https://www.gccupc.org/gupc-5201-1/30-gupc-5201-x8664-ubuntu-1204/file
# Where to get url for wget? Try git

ENV LLVM_GIT=git://github.com/Intrepid/llvm-upc.git
ENV CUPC_GIT=git://github.com/Intrepid/clang-upc.git 
#ENV CUPC_ARCHIVE=${PREFIX}/src/cupc.tar.gz
ENV CUPC_SRC_DIR=${PREFIX}/src
ENV CUPC_BLD_DIR=${PREFIX}/src/cupc/tools/bld
ENV CUPC_INSTALL_DIR=${PREFIX}
RUN mkdir -p ${CUPC_INSTALL_DIR}
RUN mkdir -p ${CUPC_SRC_DIR}
#RUN wget --quiet ${CUPC_URL} --output-document=${CUPC_ARCHIVE}
WORKDIR ${CUPC_SRC_DIR}
RUN  git clone ${LLVM_GIT} cupc >/dev/null 2>&1

WORKDIR ${CUPC_SRC_DIR}/cupc/tools
RUN  git clone ${CUPC_GIT} clang >/dev/null 2>&1
RUN mkdir -p bld
WORKDIR bld
RUN cmake  -DLLVM_TARGETS_TO_BUILD="X86" ../..
#WORKDIR ${CUPC_BLD_DIR}
#RUN ./configure --prefix=${CUPC_INSTALL_DIR} --enable-languages=c,c++
RUN make -j ${NPROCS}
#RUN make VERBOSE=1

RUN make install
RUN rm -rf ${CUPC_ARCHIVE}

ENV PATH=${PREFIX}/bin:$PATH
RUN adduser --disabled-password --gecos "" bale_user
CMD ["/bin/bash"]