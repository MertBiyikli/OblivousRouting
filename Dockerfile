# syntax=docker/dockerfile:1

########################################################
# Base image
########################################################
FROM ubuntu:24.04 AS base

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    && update-ca-certificates \
    && rm -rf /var/lib/apt/lists/*


########################################################
# Builder stage
########################################################
FROM base AS builder

# ONLY the necessary build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    ninja-build \
    git \
    pkg-config \
    wget \
    unzip \
    libeigen3-dev \
    libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt


########################################################
# Build OR-Tools v9.10 from source
########################################################
RUN git clone --depth=1 --branch v9.10 https://github.com/google/or-tools.git or-tools-src

WORKDIR /opt/or-tools-src

RUN cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
    && cmake --build build -j"$(nproc)" \
    && cmake --install build


########################################################
# Build your project
########################################################
FROM builder AS build-app

WORKDIR /app
COPY . .

RUN cmake -S . -B build -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DUSE_OPENMP=ON \
      -DBUILD_TESTS=OFF \
    && cmake --build build -j"$(nproc)"


########################################################
# Final runtime image
########################################################
FROM base AS runtime

RUN apt-get update && apt-get install -y --no-install-recommends \
    libstdc++6 \
    libgomp1 \
    libboost-program-options-dev \
    libboost-serialization-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY --from=build-app /app/build/oblivious_routing /app/oblivious_routing

ENTRYPOINT ["./oblivious_routing"]
