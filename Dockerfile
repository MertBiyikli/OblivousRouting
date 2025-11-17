# syntax=docker/dockerfile:1

FROM ubuntu:24.04 AS base
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    && update-ca-certificates \
    && rm -rf /var/lib/apt/lists/*


#########################################################
# BUILDER
#########################################################
FROM base AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential cmake ninja-build git pkg-config \
    wget unzip \
    libeigen3-dev libboost-all-dev \
    absl-dev libre2-dev zlib1g-dev \
    libprotobuf-dev protobuf-compiler \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt


#########################################################
# Build OR-Tools from source (stable & correct)
#########################################################

# Clone OR-tools (v9.10) with full history OFF
RUN git clone --depth=1 --branch v9.10 https://github.com/google/or-tools.git or-tools-src

WORKDIR /opt/or-tools-src

RUN cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
    && cmake --build build -j"$(nproc)" \
    && cmake --install build


#########################################################
# Build your project
#########################################################
FROM builder AS build-app

WORKDIR /app

COPY . .

RUN cmake -S . -B build -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DUSE_OPENMP=ON \
      -DBUILD_TESTS=OFF \
    && cmake --build build -j"$(nproc)"


#########################################################
# FINAL RUNTIME IMAGE
#########################################################
FROM base AS runtime

RUN apt-get update && apt-get install -y --no-install-recommends \
    libstdc++6 libgomp1 \
    libboost-program-options-dev libboost-serialization-dev \
    libre2-9 absl-dev zlib1g \
    libprotobuf-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY --from=build-app /app/build/oblivious_routing /app/oblivious_routing

ENTRYPOINT ["./oblivious_routing"]
