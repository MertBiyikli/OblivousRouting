# syntax=docker/dockerfile:1

FROM ubuntu:24.04 AS base

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    && update-ca-certificates \
    && rm -rf /var/lib/apt/lists/*


##############################
# Builder
##############################
FROM base AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential cmake ninja-build git pkg-config \
    wget unzip \
    libeigen3-dev libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

RUN rm -rf /opt/or-tools-src  # <-- IMPORTANT: remove cached folder

RUN git clone https://github.com/google/or-tools.git or-tools-src \
    && cd or-tools-src \
    && git checkout tags/v9.8


WORKDIR /opt/or-tools-src

RUN cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
    && cmake --build build -j"$(nproc)" \
    && cmake --install build



##############################
# Build your application
##############################
FROM builder AS build-app

WORKDIR /app
COPY . .

RUN cmake -S . -B build -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DUSE_OPENMP=ON \
      -DBUILD_TESTS=OFF \
    && cmake --build build -j"$(nproc)"


##############################
# Runtime image
##############################
FROM base AS runtime

RUN apt-get update && apt-get install -y --no-install-recommends \
    libstdc++6 libgomp1 \
    libboost-program-options-dev libboost-serialization-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY --from=build-app /app/build/oblivious_routing /app/oblivious_routing

ENTRYPOINT ["./oblivious_routing"]
