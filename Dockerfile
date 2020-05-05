FROM ubuntu:latest

RUN apt -y update && apt install -y

RUN apt -y install build-essential openmpi-bin openmpi-common

