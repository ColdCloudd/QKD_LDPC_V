FROM ubuntu:22.04 AS builder

RUN apt-get update && apt-get install -y \
    ca-certificates \
    gpg \
    wget \
    g++ \
    git \
    build-essential \
    && wget https://github.com/Kitware/CMake/releases/download/v4.1.2/cmake-4.1.2-linux-x86_64.sh \
    && mkdir /opt/cmake \
    && sh cmake-4.1.2-linux-x86_64.sh --prefix=/opt/cmake --skip-license \
    && ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake \
    && cmake --version \
    && rm -rf /var/lib/apt/lists/*\
    && rm cmake-4.1.2-linux-x86_64.sh

WORKDIR /app

COPY CMakeLists.txt /app/
COPY example/ /app/example/
COPY src/ /app/src/ 

RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
    && cmake --build build --parallel

FROM ubuntu:22.04 AS runtime

WORKDIR /app

COPY --from=builder /app/build/QKD_LDPC /app/QKD_LDPC

ENTRYPOINT ["/app/QKD_LDPC"]