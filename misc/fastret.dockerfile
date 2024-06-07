# Usage:
#   cd "$(git rev-parse --show-toplevel)/misc"
#   docker build --build-arg JDK_VERSION=8  -t "toscm/fastret:0.1.0-jdk8" -f fastret.dockerfile .
#   docker build --build-arg JDK_VERSION=11 -t "toscm/fastret:0.1.0-jdk11" -f fastret.dockerfile .
#   docker build --build-arg JDK_VERSION=17 -t "toscm/fastret:0.1.0-jdk17" -f fastret.dockerfile .
#   docker build --build-arg JDK_VERSION=18 -t "toscm/fastret:0.1.0-jdk18" -f fastret.dockerfile .
#   docker build --build-arg JDK_VERSION=19 -t "toscm/fastret:0.1.0-jdk19" -f fastret.dockerfile .
#   docker build --build-arg JDK_VERSION=21 -t "toscm/fastret:0.1.0-jdk21" -f fastret.dockerfile .
#   docker run -it --rm "toscm/fastret:0.1.0-jdk8" /bin/bash
#   docker push "toscm/fastret:0.1.0-jdk8" && docker push "toscm/fastret:0.1.0-jdk8"
FROM toscm/fastret-base:0.1.0

ARG JDK_VERSION

USER root
WORKDIR /workspace

# Copy all files individually so early steps can be cached even if later steps change.
COPY scripts/install-jdk.sh .
RUN bash install-jdk.sh "$JDK_VERSION"

COPY scripts/install-rcdk.R .
RUN Rscript install-rcdk.R

COPY scripts/test-rcdk.R .
RUN Rscript test-rcdk.R

CMD ["Rscript", "-e", "FastRet::start_gui()"]
