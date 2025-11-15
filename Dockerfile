# ============================================================
#  Stage 1: Build (compiler + dependencies + build your binary)
# ============================================================
FROM debian:stable-slim AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN echo "AMGCL_DIR=$AMGCL_DIR" \
 && echo "ORTOOLS_DIR=$ORTOOLS_DIR" \
 && ls -l $AMGCL_DIR \
 && ls -l $ORTOOLS_DIR


RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential cmake git wget curl pkg-config \
    libomp-dev libeigen3-dev unzip ca-certificates \
    && update-ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

# OR-Tools C++ SDK
WORKDIR /opt

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates && update-ca-certificates

RUN wget https://github.com/google/or-tools/releases/download/v9.10/or-tools_amd64_ubuntu-22.04_cpp_v9.10.4067.tar.gz \
    && tar -xzf or-tools_amd64_ubuntu-22.04_cpp_v9.10.4067.tar.gz \
    && extracted=$(tar -tzf or-tools_amd64_ubuntu-22.04_cpp_v9.10.4067.tar.gz | head -1 | cut -f1 -d"/") \
    && mv "$extracted" /opt/or-tools \
    && rm or-tools_amd64_ubuntu-22.04_cpp_v9.10.4067.tar.gz


ENV ORTOOLS_DIR=/opt/or-tools

# AMGCL header-only solver
RUN git clone --depth 1 https://github.com/ddemidov/amgcl.git
ENV AMGCL_DIR=/opt/amgcl

# Copy project AFTER dependencies (max caching)
WORKDIR /build
COPY . .

RUN cmake -S . -B /build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$ORTOOLS_DIR" \
    -DCMAKE_CXX_FLAGS="-O3 -march=native -fopenmp" \
    2>&1 | tee /tmp/cmake_config.log \
 && cmake --build /build -j$(nproc) \
    2>&1 | tee /tmp/cmake_build.log \
 && echo "==== BUILD DIRECTORY CONTENT ====" \
 && ls -R /build






# ============================================================
#  Stage 2: Runtime (small)
# ============================================================
FROM debian:stable-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential cmake git wget curl pkg-config \
    libomp-dev libeigen3-dev unzip ca-certificates \
    lib


COPY --from=builder /opt/or-tools /opt/or-tools
ENV ORTOOLS_DIR=/opt/or-tools

COPY --from=builder /build/oblivious_routing /usr/local/bin/oblivious_routing


WORKDIR /data
ENTRYPOINT ["oblivious_routing"]
