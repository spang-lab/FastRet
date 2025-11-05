# Usage:
#   $env:DOCKER_BUILDKIT=1
#   cd "$(git rev-parse --show-toplevel)"
#   docker build -t "toscm/fastret-base:1.2.2" -t "toscm/fastret-base:latest" -f misc/fastret-base.dockerfile .
#   docker run -it --rm -p 8080:8080 "toscm/fastret-base:1.2.2"
#   docker run -it --rm -p 8080:8080 -v "$(pwd)/..:/workspace/FastRet" "toscm/fastret-base:1.2.2" /bin/bash
#   docker push "toscm/fastret-base:1.2.2" && docker push "toscm/fastret-base:latest"
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Berlin

USER root
WORKDIR /workspace/FastRet/misc/scripts

# Copy all files individually so early steps can be cached even if later steps change.
COPY ./misc/scripts/install-sys-pkgs.sh /workspace/FastRet/misc/scripts
RUN bash install-sys-pkgs.sh

COPY ./misc/scripts/install-r.sh /workspace/FastRet/misc/scripts
RUN bash install-r.sh

COPY ./misc/scripts/install-r-dev-pkgs.R /workspace/FastRet/misc/scripts
RUN Rscript install-r-dev-pkgs.R

RUN R CMD javareconf

COPY ./misc/scripts/install-rjava.R /workspace/FastRet/misc/scripts
RUN Rscript install-rjava.R

COPY ./misc/scripts/install-radian.sh /workspace/FastRet/misc/scripts
RUN bash install-radian.sh

COPY ./misc/scripts/install-en-us-locales.sh /workspace/FastRet/misc/scripts
RUN bash install-en-us-locales.sh

COPY ./misc/scripts/configure-users.sh /workspace/FastRet/misc/scripts
RUN bash configure-users.sh

COPY . /workspace/FastRet
RUN Rscript install-fastret.R --local --verbose

USER shiny
WORKDIR /home/shiny

CMD ["Rscript", "-e", "FastRet::start_gui()"]
