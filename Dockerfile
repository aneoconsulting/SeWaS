FROM ubuntu:18.04

RUN DEBIAN_FRONTEND=noninteractive apt-get -y update
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y

RUN DEBIAN_FRONTEND=noninteractive apt-get -y install build-essential openmpi-bin openmpi-common

COPY public/thirdparty/install /usr
COPY public/sewas /home/sewas