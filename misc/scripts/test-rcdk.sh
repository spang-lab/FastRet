#!/usr/bin/bash

# Usage: bash test-rcdk.sh [--verbose]

echo "------------------------------------------"
for jdk in default-jdk openjdk-8-jdk openjdk-11-jdk openjdk-17-jdk openjdk-18-jdk openjdk-19-jdk openjdk-21-jdk
do
    echo "Testing $jdk ..."
    bash ./install-jdk.sh "$jdk" $1
    Rscript ./install-rcdk.R $$1
    Rscript ./test-rcdk.R $1
    echo "------------------------------------------"
done
