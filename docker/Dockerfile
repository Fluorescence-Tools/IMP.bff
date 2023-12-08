FROM ubuntu:22.04
LABEL authors="tpeulen"

RUN apt update && DEBIAN_FRONTEND=noninteractive apt install -y tzdata
RUN apt install -y vim build-essential git cmake net-tools gdb clang

RUN DEBIAN_FRONTEND="noninteractive" apt-get update \
    && apt-get install -y ninja-build \
    && apt-get clean

WORKDIR /work
#ENTRYPOINT ["top", "-b"]