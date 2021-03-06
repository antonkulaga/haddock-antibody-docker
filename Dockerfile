FROM mambaorg/micromamba:latest
#non-root user of the mambaorg/micromamba is buggy, so we fallback to root
USER root
RUN apt-get update && apt-get install -y git wget acl libxml2-dev build-essential libcurl4-openssl-dev pkg-config m4 libtool automake autoconf libjson-c-dev libgl1-mesa-dev libgl-dev libc++-11-dev libc++abi-11-dev
ENV BASE=/opt
WORKDIR $BASE
RUN git clone https://github.com/haddocking/haddock-tools
WORKDIR $BASE/haddock-tools
RUN make
WORKDIR $BASE
RUN git clone https://github.com/antonkulaga/HADDOCK-antibody-antigen.git haddock-antibody
ENV PATH="$PATH:$BASE/haddock-antibody"
COPY environment.yaml $BASE/environment.yaml
RUN micromamba env create -y -f $BASE/environment.yaml
ENV ENV_NAME=haddock-antibody
ENV PATH="$PATH:$BASE"
RUN mkdir /data && chmod -R 777 /data
RUN chown root:root /usr/local/bin/_entrypoint.sh && chmod +x /usr/local/bin/_entrypoint.sh
VOLUME ["/data"]
WORKDIR /data
COPY extract_fasta.py $BASE/extract_fasta.py
COPY start.py $BASE/start.py
CMD python $BASE/start.py
