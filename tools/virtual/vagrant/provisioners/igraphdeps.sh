#! /bin/bash

apt-get -y install git autoconf automake bison flex libtool \
    libxml2-dev libgmp-dev docbook2x source-highlight libxml2-utils \
    r-base-core
apt-get clean
