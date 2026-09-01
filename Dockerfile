FROM ubuntu:18.04

RUN DEBIAN_FRONTEND=noninteractive apt-get -y update
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y
RUN DEBIAN_FRONTEND=noninteractive apt-get -y install build-essential libboost-program-options-dev openmpi-bin openmpi-common

# Install thirdparty libraries
COPY ./public/install /usr

# Install sewas artifacts
COPY ./public/sewas /home/sewas
RUN chmod u+rwx /home/sewas/sewas
