#!/usr/bin/bash
echo "------------------------------------------"
for jdk in default-jdk openjdk-8-jdk openjdk-11-jdk openjdk-17-jdk openjdk-18-jdk openjdk-19-jdk openjdk-21-jdk
do
    echo "Testing $jdk ..."
    bash ./install_jdk.sh "$jdk"
    Rscript ./install_rcdk.R
    Rscript ./test_rcdk.R
    echo "------------------------------------------"
done
