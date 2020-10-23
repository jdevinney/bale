
FROM ubuntu:18.04

MAINTAINER bale <bale@super.org>

# Build-arg: NPROCS for parallel make
ARG NPROCS=8

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
	flex \
	gcc \
	g++ \
	gfortran \
	gcc-multilib \
        git \
	libnuma-dev \
        libtool \
        m4 \
        make \
	openssh-client \
        patch \
	pkg-config \
	wget

ENV OFI_GIT=https://github.com/ofiwg/libfabric
#RUN mkdir -p ${PREFIX}/OFI
RUN git clone ${OFI_GIT} ${PREFIX}/OFI  >/dev/null 2>&1
WORKDIR ${PREFIX}/OFI
RUN ./autogen.sh
RUN ./configure
RUN make -j${NPROCS}
RUN make install

WORKDIR ${PREFIX}
ENV HYDRA_URL=http://www.mpich.org/static/downloads/3.3.2/hydra-3.3.2.tar.gz
ENV HYDRA_ARCHIVE=${PREFIX}/hydra-3.3.2.tar.gz
ENV HYDRA_VERSION=hydra-3.3.2
RUN wget --quiet ${HYDRA_URL} 
RUN tar -zxf ${HYDRA_ARCHIVE} 
WORKDIR ${PREFIX}/hydra-3.3.2
#RUN ./autogen.sh
RUN ./configure
RUN make -j${NPROCS}
RUN make install
RUN rm -f $HYDRA_ARCHIVE}

WORKDIR ${PREFIX}
ENV SOS_GIT=https://github.com/Sandia-OpenSHMEM/SOS
RUN git clone ${SOS_GIT} ${PREFIX}/SOS  >/dev/null 2>&1
WORKDIR ${PREFIX}/SOS
RUN ./autogen.sh
RUN ./configure --with-ofi=${PREFIX} --disable-fortran --enable-pmi-simple --disable-error-checking --enable-memcpy
RUN make -j${NPROCS}
RUN make install

WORKDIR ${PREFIX}
RUN export PATH=${PREFIX}/bin:$PATH
RUN export LD_LIBRARY_PATH=${PREFIX}/lib:$LD_LIBRARY_PATH
RUN adduser --disabled-password --gecos "" bale_user
CMD ["/bin/bash"]
