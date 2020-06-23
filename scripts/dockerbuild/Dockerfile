# creates an image containing vg and hal2vg

# build on compatible vg image
FROM quay.io/vgteam/vg:v1.24.0

# update system and install dependencies not present in vg image
RUN apt-get -qq update && apt-get -qq install -y libhdf5-dev build-essential python3-dev python3-pip git cmake

# copy current directory to docker

RUN git clone https://github.com/glennhickey/hal2vg.git --recursive && cd /hal2vg && make

# add hal2vg to the PATH
ENV PATH /hal2vg:/hal2vg/deps/hal/bin:$PATH
