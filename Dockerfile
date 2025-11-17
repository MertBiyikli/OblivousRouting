# ============================================================
# Stage 1: Build OR-Tools + your binary
# ============================================================
FROM ubuntu:24.04 AS builder

ARG DEBIAN_FRONTEND=noninteractive

# Fix SSL trust for git/github
RUN apt-get update && apt-get install -y ca-certificates && update-ca-certificates

# System build deps
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential cmake ninja-build git pkg-config \
    wget unzip \
    libeigen3-dev libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

# ------------------------------------------------------------
# 1. Download & build OR-Tools (C++)
# ------------------------------------------------------------
ENV ORTOOLS_VERSION=v9.10

ARG ORTOOLS_SOURCE=/opt/or-tools

# ----------------------------------------
# OR-Tools caching support
# ----------------------------------------

# Create cache directory so COPY never fails
RUN mkdir -p /opt/or-tools-source

# Copy any cached OR-Tools source (may be empty)
COPY .cache/ortools/ /opt/or-tools-source/

# Build OR-Tools from cached source if available, otherwise download
RUN if [ -d "/opt/or-tools-source" ] && [ "$(ls -A /opt/or-tools-source)" ]; then \
      echo "Using cached OR-Tools source"; \
      cd /opt/or-tools-source && cmake -S . -B build -G Ninja && cmake --build build -j"$(nproc)" && cmake --install build; \
    else \
      echo "Downloading OR-Tools from GitHub"; \
      cd /opt && git clone --depth=1 https://github.com/google/or-tools.git or-tools-source && \
      cd or-tools-source && cmake -S . -B build -G Ninja && cmake --build build -j"$(nproc)" && cmake --install build; \
    fi


WORKDIR /opt/or-tools

RUN cmake -S . -B build -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DBUILD_DEPS=ON \
      -DBUILD_CXX=ON \
      -DBUILD_PYTHON=OFF \
      -DBUILD_JAVA=OFF \
      -DBUILD_EXAMPLES=OFF \
      -DBUILD_TESTING=OFF \
 && cmake --build build -j"$(nproc)" \
 && cmake --install build

# ------------------------------------------------------------
# 2. Build your project (NO TESTS INSIDE DOCKER)
# ------------------------------------------------------------
WORKDIR /app
COPY . .

RUN cmake -S . -B build -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DUSE_OPENMP=ON \
      -DBUILD_TESTS=OFF \
 && cmake --build build -j"$(nproc)"

# ============================================================
# Stage 2: Runtime image
# ============================================================
FROM ubuntu:24.04 AS runtime

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    libstdc++6 libgomp1 libboost-program-options-dev libboost-serialization-dev \
    && rm -rf /var/lib/apt/lists/*

# App binary
COPY --from=builder /app/build/oblivious_routing /usr/local/bin/oblivious_routing

WORKDIR /data

ENTRYPOINT ["oblivious_routing"]
CMD ["--help"]
