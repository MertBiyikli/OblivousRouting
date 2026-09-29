# Reproducible application build.
#   docker build --build-arg TOOLCHAIN_IMAGE=ghcr.io/...@sha256:... .
ARG TOOLCHAIN_IMAGE=ghcr.io/mertbiyikli/ortools@sha256:9f9e38fa700e5bdea953bd1ca43efe5d0670e4a153d877b659a2e83ada6495d7

FROM ${TOOLCHAIN_IMAGE} AS builder
WORKDIR /src
COPY . .

RUN cmake --preset release \
    && cmake --build --preset release

# For now use the same pinned toolchain environment for runtime correctness.
# A separate slim runtime image can be added later once runtime shared-library
# dependencies are explicitly packaged and tested.
FROM ${TOOLCHAIN_IMAGE} AS runtime
WORKDIR /app
COPY --from=builder /src/build/release/oblivious_routing /usr/local/bin/oblivious_routing

ENTRYPOINT ["/usr/local/bin/oblivious_routing"]
