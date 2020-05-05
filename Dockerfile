FROM ubuntu:18.04

RUN DEBIAN_FRONTEND=noninteractive apt-get -y update
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y

RUN DEBIAN_FRONTEND=noninteractive apt-get -y install build-essential libboost-program-options-dev openmpi-bin openmpi-common

COPY ./public/sewas /home/sewas
COPY ./public/install /usr
