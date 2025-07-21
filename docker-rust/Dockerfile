FROM rust:1.88

RUN apt-get update
RUN apt-get install -y gcc-x86-64-linux-gnu g++-x86-64-linux-gnu clang libclang-dev
RUN rustup component add clippy
RUN rustup update
RUN cargo install cargo-license
RUN git clone http://bitbucket.pacificbiosciences.com:7990/scm/~mholt/license_checker.git


# set up python 3.10 via micromamba to enable python + pip in SVTopoVz testing
RUN wget https://micromamba.snakepit.net/api/micromamba/linux-64/latest -O micromamba.tar.bz2
RUN mkdir -p /opt/micromamba
RUN tar -xvjf micromamba.tar.bz2 -C /opt/micromamba --strip-components=1
RUN ln -s /opt/micromamba/micromamba /usr/local/bin/micromamba
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
RUN micromamba create -y -n py310 "python=3.10" "pytest=8.3.5"
