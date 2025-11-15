# ============================================================
# Stage 1: Build OR-Tools + your binary
# ============================================================
FROM ubuntu:24.04 AS builder

ARG DEBIAN_FRONTEND=noninteractive

# --- Fix TLS/SSL issues inside Docker on macOS/Ubuntu ---
RUN apt-get update && apt-get install -y ca-certificates && update-ca-certificates


RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential cmake ninja-build git pkg-config \
    wget unzip \
    libeigen3-dev libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

# ------------------------------------------------------------
# 1. Download & build OR-Tools C++ (stable version)
# ------------------------------------------------------------
ENV ORTOOLS_VERSION v9.10

RUN git clone --branch ${ORTOOLS_VERSION} --depth=1 https://github.com/google/or-tools.git

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
# 2. Build your project
# ------------------------------------------------------------
WORKDIR /app
COPY . .

RUN cmake -S . -B build -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DUSE_OPENMP=ON \
 && cmake --build build -j"$(nproc)"

RUN cmake --build build -j"$(nproc)" \
 && cd build && ctest --show-only


# ============================================================
# Stage 2: Runtime image
# ============================================================
FROM ubuntu:24.04 AS runtime

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    libstdc++6 libgomp1 libboost-program-options-dev libboost-serialization-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/local/bin
COPY --from=builder /usr/local/bin/ /usr/local/bin/
COPY --from=builder /app/build/oblivious_routing /usr/local/bin/oblivious_routing


WORKDIR /app/build
COPY --from=builder /app/build/ /app/build/


WORKDIR /data

ENTRYPOINT ["oblivious_routing"]
CMD ["--help"]
