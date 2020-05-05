FROM ubuntu:latest

RUN DEBIAN_FRONTEND=noninteractive apt-get -y update
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y

RUN DEBIAN_FRONTEND=noninteractive apt-get -y install build-essential openmpi-bin openmpi-common

COPY . /home/SeWaS