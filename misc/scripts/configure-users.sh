#!/bin/bash

useradd --create-home shiny
usermod -a -G staff shiny # add shiny user to staff group so they can install into /usr/local/lib/R/site-library
chsh -s /bin/bash shiny
