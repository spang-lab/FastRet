#!/bin/bash

# Required to prevent the following R CMD check warning:
# >>> checking R files for syntax errors ...
# >>> Warning in Sys.setlocale("LC_CTYPE", "en_US.UTF-8") :
# >>>   OS reports request to set locale to "en_US.UTF-8" cannot be honored

# Check if the script is run as root
if [ "$EUID" -ne 0 ]
  then echo "Please run as root"
  exit
fi

# Generate en_US.UTF-8 locale
apt-get update
apt-get install --no-install-recommends -y locales
apt-get clean
rm -rf /var/lib/apt/lists/*
sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen
locale-gen
