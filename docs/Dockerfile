# Base image with Linux (Debian)
FROM ubuntu:22.04

# Set non-interactive mode for apt-get
ENV DEBIAN_FRONTEND=noninteractive

# Install GCC 10, R, Python, and necessary dependencies
RUN apt-get update \
    && apt-get install -y \
        build-essential \
        cmake \
        curl \
        g++-10 \
        gcc-10 \
        gfortran-10 \
        libatlas-base-dev \
        libc6 \
        libcurl4-openssl-dev \
        libhdf5-dev \
        liblapack-dev \
        libopenblas-dev \
        libssl-dev \
        python3 \
        python3-dev \
        python3-pip \
        r-base \
        software-properties-common \
        wget \
        zlib1g-dev \
        unzip \
        bash \
    && apt-get install -y unzip bash \
    && ln -sf /bin/bash /bin/sh \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set GCC 10 as default (update-alternatives to switch versions if needed)
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10 100 \
    && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-10 100 \
    && update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-10 100
ENV CC=gcc-10
ENV CXX=g++-10

# Verify installations
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN gcc --version && python3 --version && R --version

# Install bgenix and compile
RUN wget http://code.enkre.net/bgen/tarball/release/bgen.tgz \
        -O /usr/bin/bgen.tgz \
    && tar -xzf /usr/bin/bgen.tgz \
        -C /usr/bin \
    && mv /usr/bin/bgen.tgz /usr/bin/bgen_dir \
    && cd /usr/bin/bgen_dir \
    && CXX=g++-10 CC=gcc-10 ./waf configure \
    && ./waf \
    && ln -s /usr/bin/bgen_dir/build/apps/bgenix /usr/bin/bgenix
    
# Install regenie and compile 
RUN wget https://github.com/rgcgithub/regenie/archive/refs/heads/master.zip \
        -O /usr/bin/regenie.zip \
    && unzip /usr/bin/regenie.zip \
        -d /usr/bin/regenie_dir \
    && rm /usr/bin/regenie.zip \
    && mkdir /usr/bin/regenie_dir/regenie-master/build \
    && cd /usr/bin/regenie_dir/regenie-master/build \
    && BGEN_PATH=/usr/bin/bgen_dir cmake .. \
    && make \
    && ln -s /usr/bin/regenie_dir/regenie-master/build/regenie /usr/bin/regenie

# Install plink2
RUN wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20241203.zip \
        -O /usr/bin/plink2.zip \
    && unzip /usr/bin/plink2.zip \
        -d /usr/bin/plink2_dir \
    && rm /usr/bin/plink2.zip \
    && ln -s /usr/bin/plink2_dir/plink2 /usr/bin/plink2

# Add R packages 
RUN R -e "install.packages(c('optparse', 'RColorBrewer'), repos='https://cloud.r-project.org')"
RUN PIP_PATH=$(pip show agentp | grep Location | awk '{print $2"/agentp"}') \
    && R -e "install.packages('$PIP_PATH/external/plink2R', repos=NULL, type='source')"

# Set working directory in the container
WORKDIR /usr/home

# Download example data
RUN wget https://vanderbilt.box.com/shared/static/o5s5bwwey0clt8ktabqg1abmubg2t7in.gz \
        -O /usr/home/example.tar.gz \
    && tar -xzf /usr/home/example.tar.gz \
    && rm /usr/home/example.tar.gz

# Install agentp 
RUN pip3 uninstall agentp
RUN pip3 install agentp

# Download TWAS models to agentp path 
#RUN PIP_PATH=$(pip show agentp | grep Location | awk '{print $2"/agentp"}') \
#    && wget https://vanderbilt.box.com/shared/static/fs4rltegomfq568y6vyjxq22vizzv48u.gz \
RUN PIP_PATH=$(pip show agentp | grep Location | awk '{print $2"/agentp"}') \
    && wget https://vanderbilt.box.com/shared/static/10hz24rk7z9r6oh7h3st84vicv6ksqfq.gz \
        -O "$PIP_PATH/models/models.tar.gz" \
    && tar -xzf "$PIP_PATH/models/models.tar.gz" \
        --strip-components=1 \
        -C "$PIP_PATH/models" \
    && rm "$PIP_PATH/models/models.tar.gz"

