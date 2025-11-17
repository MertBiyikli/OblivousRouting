# syntax=docker/dockerfile:1

###############################################
# Base image (shared between build + runtime)
###############################################
FROM ubuntu:24.04 AS base

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    && update-ca-certificates \
    && rm -rf /var/lib/apt/lists/*


###############################################
# Builder image
###############################################
FROM base AS builder

# Build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential cmake ninja-build git pkg-config \
    wget unzip \
    libeigen3-dev libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt


###########################################################
# OR-Tools installation using PRECOMPILED binary + cache
###########################################################

# Ensure the target directory exists to avoid COPY errors
RUN mkdir -p /opt/or-tools

# Bring in cached precompiled OR-Tools from GitHub Actions
COPY .cache/ortools/ /opt/or-tools/

# If cache is empty, download OR-Tools for Ubuntu 24.04
RUN if [ ! -f "/opt/or-tools/lib/libortools.a" ]; then \
      echo "Downloading precompiled OR-Tools v9.10 (generic Linux x86_64)..."; \
      wget -q https://github.com/google/or-tools/releases/download/v9.10/or-tools_amd64_linux_v9.10.4067.tar.gz -O or-tools.tar.gz && \
      mkdir -p /opt/or-tools && \
      tar -xzf or-tools.tar.gz --strip-components=1 -C /opt/or-tools && \
      rm or-tools.tar.gz; \
    else \
      echo "Using cached OR-Tools."; \
    fi


# Make OR-Tools available via CMake
ENV CMAKE_PREFIX_PATH=/opt/or-tools
ENV LD_LIBRARY_PATH=/opt/or-tools/lib:${LD_LIBRARY_PATH}


###############################################
# Build your project
###############################################
WORKDIR /app

# Copy entire repository including extern/
COPY . .

RUN cmake -S . -B build -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DUSE_OPENMP=ON \
      -DBUILD_TESTS=OFF \
    && cmake --build build -j"$(nproc)"


###############################################
# Runtime image (small and clean)
###############################################
FROM base AS runtime

RUN apt-get update && apt-get install -y --no-install-recommends \
    libstdc++6 libgomp1 \
    libboost-program-options-dev libboost-serialization-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy binary from build stage
COPY --from=builder /app/build/oblivious_routing /app/oblivious_routing

ENTRYPOINT ["./oblivious_routing"]
