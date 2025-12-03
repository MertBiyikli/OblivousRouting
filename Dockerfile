FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# System packages needed for OR-Tools
RUN apt-get update && apt-get install -y \
    git cmake build-essential ninja-build \
    libz-dev libbz2-dev libeigen3-dev \
    libgmp-dev libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*

# Build OR-Tools from source (stable v9.6)
RUN git clone https://github.com/google/or-tools.git --branch v9.6 --depth 1
RUN mkdir ortools-build && cd ortools-build && \
    cmake ../or-tools -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DBUILD_DEPS=ON \
      -DBUILD_SAMPLES=OFF \
    && cmake --build . -j"$(nproc)" \
    && cmake --install . --prefix /usr/local

# Test installation
RUN ls -R /usr/local/lib/cmake/ortools

# Final image has OR-Tools only
