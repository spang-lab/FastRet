# Usage:
#   cd "$(git rev-parse --show-toplevel)/misc"
#   docker build -t toscm/rcdk-test:0.1.0 -f rcdk-testttt.dockerfile .
#   docker run -it --rm toscm/rcdk-test:0.1.0
#   docker run -it --rm toscm/rcdk-test:0.1.0 /bin/bash
#   docker run -it --rm toscm/rcdk-test:0.1.0 /bin/bash test-rcdk.sh --verbose
#   docker push toscm/rcdk-test:0.1.0
FROM ubuntu:22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Berlin

# Install R as described at https://cran.r-project.org/bin/linux/ubuntu/index.html
RUN apt-get update
RUN apt-get install --no-install-recommends -y software-properties-common dirmngr wget
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc >> /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN apt-get install --no-install-recommends -y gpg-agent
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
RUN add-apt-repository "ppa:c2d4u.team/c2d4u4.0+"
RUN apt-get update
RUN apt-get install --no-install-recommends -y libbz2-dev libdeflate-dev r-base r-base-dev

COPY scripts /workspace
WORKDIR /workspace

CMD ["/bin/bash", "test-rcdk.sh"]
