FROM ubuntu:latest

RUN apt-get -y update && apt-get install -y

RUN apt-get -y install build-essential openmpi-bin openmpi-common

