##############################
# Stage 1 — Build OR-Tools
##############################
FROM ubuntu:24.04 AS ortools-builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    git cmake ninja-build build-essential \
    libz-dev libbz2-dev libeigen3-dev \
    libgmp-dev libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*

# Build OR-Tools v9.6
RUN git clone https://github.com/google/or-tools.git --branch v9.6 --depth 1
RUN mkdir ortools-build && cd ortools-build && \
    cmake ../or-tools -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DBUILD_SAMPLES=OFF \
      -DBUILD_DEPS=ON \
    && cmake --build . -j"$(nproc)" \
    && cmake --install . --prefix /opt/ortools


##############################
# Stage 2 — Build your project
##############################
FROM ubuntu:24.04 AS cpp-builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    cmake ninja-build build-essential \
    libz-dev libbz2-dev libeigen3-dev \
    && rm -rf /var/lib/apt/lists/*

COPY --from=ortools-builder /opt/ortools /opt/ortools
ENV CMAKE_PREFIX_PATH=/opt/ortools

WORKDIR /src
COPY . .

RUN cmake -S . -B build -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DUSE_OPENMP=ON \
  && cmake --build build -j"$(nproc)"


##############################
# Stage 3 — Final tiny runtime
##############################
FROM ubuntu:24.04 AS runtime

RUN apt-get update && apt-get install -y \
    libz-dev libbz2-dev libeigen3-dev \
    && rm -rf /var/lib/apt/lists/*

COPY --from=cpp-builder /src/build/oblivious_routing /usr/local/bin/oblivious_routing
COPY --from=ortools-builder /opt/ortools /opt/ortools

ENV LD_LIBRARY_PATH=/opt/ortools/lib

ENTRYPOINT ["/usr/local/bin/oblivious_routing"]
