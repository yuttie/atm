FROM gcc:10
RUN apt-get update && apt-get install -y \
    cmake \
    libboost-chrono1.67-dev \
    libboost-timer1.67-dev \
    libboost-system1.67-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /usr/src/myapp
