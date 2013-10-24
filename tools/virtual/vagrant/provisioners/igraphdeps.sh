#! /bin/bash

apt-get -y install git autoconf automake bison flex libtool \
    libxml2-dev libgmp-dev docbook2x source-highlight libxml2-utils \
    r-base-core tcl8.5-dev tk8.5-dev \
    mesa-common-dev libglu1-mesa-dev texlive-base \
    texlive-latex-recommended texlive-fonts-extra texlive-latex-extra \
    texlive-fonts-recommended
apt-get clean
