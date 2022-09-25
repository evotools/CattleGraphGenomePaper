FROM continuumio/miniconda3 AS build

LABEL authors="andrea.talenti@ed.ac.uk" \
      description="Docker image containing base requirements for GraphSeq pipelines"

RUN apt-get update && \
  apt-get install -y gcc g++ make cmake doxygen && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install the package as normal:
COPY environment.yml .

# Install mamba to speed up the process
RUN conda install -c conda-forge -y mamba

# Create the environment
RUN mamba env create -f environment.yml

# Install conda-pack:
RUN mamba install -c conda-forge conda-pack

# Use conda-pack to create a standalone enviornment
# in /venv:
RUN conda-pack -n graphseq -o /tmp/env.tar && \
  mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
  rm /tmp/env.tar

# We've put venv in same path it'll be in final image,
# so now fix up paths:
RUN /venv/bin/conda-unpack


# The runtime-stage image; we can use Debian as the
# base image since the Conda env also includes Python
# for us.
FROM debian:buster AS runtime

# Install procps in debian to make it compatible with reporting
RUN apt-get update && \
  apt install -y procps file wget && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Copy /venv from the previous stage:
COPY --from=build /venv /venv

# When image is run, run the code with the environment
# activated:
ENV PATH /venv/bin/:$PATH

# prepare entry point
SHELL ["/bin/bash", "-c"]
