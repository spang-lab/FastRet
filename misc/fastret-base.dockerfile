# Usage:
#   $env:DOCKER_BUILDKIT=1
#   cd "$(git rev-parse --show-toplevel)/misc"
#   docker build -t "toscm/fastret-base:0.1.0" -t "toscm/fastret-base:latest" -f fastret-base.dockerfile .
#   docker run -it --rm -p 8080:8080 "toscm/fastret-base:latest"
#   docker run -it --rm -p 8080:8080 -v "$(pwd)/..:/workspace/FastRet" "toscm/fastret-base:latest" /bin/bash
#   docker push "toscm/fastret-base:0.1.0" && docker push "toscm/fastret-base:latest"
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Berlin

USER root
WORKDIR /workspace

# Copy all files individually so early steps can be cached even if later steps change.
COPY scripts/install-sys-pkgs.sh .
RUN bash install-sys-pkgs.sh

COPY scripts/install-r.sh .
RUN bash install-r.sh

COPY scripts/install-r-dev-pkgs.R .
RUN Rscript install-r-dev-pkgs.R

RUN R CMD javareconf

COPY scripts/install-rjava.R .
RUN Rscript install-rjava.R

COPY scripts/install-fastret.R .
RUN Rscript install-fastret.R --branch fix-actions --verbose

CMD ["Rscript", "-e", "FastRet::start_gui()"]
