############################
# Stage 1: Builder (amd64)
############################
FROM --platform=linux/amd64 ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    wget \
    ca-certificates \
    libstdc++6 \
    libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*

# Install OR-Tools (amd64)
RUN wget -q https://github.com/google/or-tools/releases/download/v9.9/or-tools_amd64_ubuntu-22.04_cpp_v9.9.3963.tar.gz && \
    tar -xzf or-tools_amd64_ubuntu-22.04_cpp_v9.9.3963.tar.gz && \
    cp -r or-tools_*/* /usr/local/ && \
    ldconfig && \
    rm -rf or-tools_*

WORKDIR /build
COPY . .

RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_OPENMP=OFF && \
    cmake --build build -j$(nproc)

############################
# Stage 2: Runtime (amd64)
############################
FROM --platform=linux/amd64 ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    libstdc++6 \
    libgomp1 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Copy OR-Tools runtime libs
COPY --from=builder /usr/local/lib /usr/local/lib
RUN ldconfig

WORKDIR /app
COPY --from=builder /build/build/oblivious_routing .

ENTRYPOINT ["./oblivious_routing"]
