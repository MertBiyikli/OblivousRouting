# syntax=docker/dockerfile:1

# -------------------------
# Base image (shared)
# -------------------------
FROM ubuntu:24.04 AS base

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    && update-ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# -------------------------
# Builder image
# -------------------------
FROM base AS builder

# Build tools + libs needed for building your project and OR-Tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential cmake ninja-build git pkg-config \
    wget unzip \
    libeigen3-dev libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

# ----------------------------------------
# OR-Tools: use cached source if provided
# ----------------------------------------

# Ensure target directory exists so COPY never errors
RUN mkdir -p /opt/or-tools-source

# Copy cached OR-Tools source from build context (from .cache/ortools)
# This directory is prepared & cached by GitHub Actions.
COPY .cache/ortools/ /opt/or-tools-source/

# If cache is non-empty, build from there; otherwise clone fresh.
RUN if [ -d "/opt/or-tools-source" ] && [ "$(ls -A /opt/or-tools-source)" ]; then \
      echo "Using cached OR-Tools source in /opt/or-tools-source"; \
    else \
      echo "No OR-Tools source in context, cloning from GitHub..."; \
      rm -rf /opt/or-tools-source && \
      git clone --depth=1 --branch v9.10 https://github.com/google/or-tools.git /opt/or-tools-source; \
    fi && \
    cd /opt/or-tools-source && \
    cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release && \
    cmake --build build -j"$(nproc)" && \
    cmake --install build

# -------------------------
# Build your project
# -------------------------
WORKDIR /app

# Copy your entire repo (including extern/, CMakeLists.txt, src/, etc.)
COPY . .

# Configure & build (tests disabled to keep image smaller / faster)
RUN cmake -S . -B build -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DUSE_OPENMP=ON \
      -DBUILD_TESTS=OFF \
    && cmake --build build -j"$(nproc)"

# -------------------------
# Runtime image
# -------------------------
FROM base AS runtime

# Only runtime dependencies (smaller image)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libstdc++6 libgomp1 \
    libboost-program-options-dev libboost-serialization-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy the built binary from the builder stage
COPY --from=builder /app/build/oblivious_routing /app/oblivious_routing

ENTRYPOINT ["./oblivious_routing"]
