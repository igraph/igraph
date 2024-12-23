FROM ubuntu:18.04
# FROM ubuntu:latest

RUN apt-get update && apt-get install -y \
    build-essential \
    python3-pip \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

COPY . /infomap/

WORKDIR /infomap

RUN make -j

RUN pip3 --no-cache-dir install --index-url https://test.pypi.org/simple/ infomap

ENTRYPOINT ["/infomap/Infomap"]
CMD ["--help"]