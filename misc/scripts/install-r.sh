#!/bin/bash

# Check if the script is run as root
if [ "$EUID" -ne 0 ]
  then echo "Please run as root"
  exit
fi

# Install R
apt-get update
apt-get install --no-install-recommends -y software-properties-common dirmngr wget
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc >> /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
apt-get install --no-install-recommends -y gpg-agent
add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
add-apt-repository "ppa:c2d4u.team/c2d4u4.0+"
apt-get update
apt-get install --no-install-recommends -y r-base r-base-dev
