FROM rocker/rstudio:3

RUN apt-get update && apt-get install -y -q \
    build-essential \
    libpcre3-dev \
    autoconf \
    automake \
    libtool \
    bison \
    git \
    libboost-dev \
    python3-dev \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/swig/swig/archive/rel-4.0.1.tar.gz \
&& tar -zxf rel-4.0.1.tar.gz \
&& cd swig-rel-4.0.1 \
&& rm -f ../rel-4.0.1.tar.gz \
&& ./autogen.sh \
&& ./configure \
&& make \
&& make install

COPY . /home/rstudio/

WORKDIR /home/rstudio

RUN make R

WORKDIR /home/rstudio/build/R/Infomap

# ENTRYPOINT ["/home/rstudio/Infomap"]
# CMD ["--help"]
