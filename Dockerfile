# Build base
FROM condaforge/mambaforge:latest AS build-base
RUN apt-get update && apt-get install -y \
      bash \
      build-essential


FROM build-base as build-conda
COPY ./environment.yml /environment.yml
RUN --mount=type=cache,target=/opt/conda/pkgs mamba env create -n babachi --file /environment.yml && echo 'conda activate babachi' >> ~/.bashrc
SHELL ["conda", "run", "--no-capture-output", "-n", "babachi", "/bin/bash", "-c"]

RUN git clone https://github.com/vierstralab/WASP && cd WASP/snp2h5 && make

#######################
# Final image
FROM condaforge/mambaforge:latest AS wasp-realigning
RUN apt-get update && apt-get install -y bash


COPY --from=build-conda /opt/conda /opt/conda
COPY --from=build-conda /WASP/ /opt/WASP
ENV PATH=/opt/conda/envs/babachi/bin:$PATH
ARG PATH=/opt/conda/envs/babachi/bin:$PATH
