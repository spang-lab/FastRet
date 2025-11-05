#!/bin/bash

# Check if the script is run as root
if [ "$EUID" -ne 0 ]
  then echo "Please run as root"
  exit
fi

# Install system packages
apt-get update
apt-get install --no-install-recommends -y curl default-jdk devscripts gdebi-core ghostscript git htop libbz2-dev libcairo2-dev libcurl4-openssl-dev libdeflate-dev libfontconfig1-dev libfreetype6-dev libfribidi-dev libgit2-dev libharfbuzz-dev libicu-dev libjpeg-dev libpng-dev libssl-dev libtiff-dev libwebp-dev libxml2-dev make micro nginx openssh-server pandoc python3-pip qpdf rsyslog supervisor texlive texlive-fonts-extra tmux wget zlib1g-dev
apt-get clean
rm -rf /var/lib/apt/lists/*
