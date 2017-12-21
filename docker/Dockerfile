# A base Docker image that other Chaste images can build from.
#
# Handles the base Ubuntu install, installing all the Chaste dependencies,
# and creating a default chaste user.
#
# Several ideas taken from the FEniCS project docker setup
# (https://bitbucket.org/fenics-project/docker).

#FROM ubuntu:trusty
FROM phusion/baseimage:latest
MAINTAINER Chaste Developers <chaste-admin@maillist.ox.ac.uk>

# Install the Chaste repo list, and dependencies metapackage
#USER root
COPY chaste.list /etc/apt/sources.list.d
RUN apt-key adv --recv-keys --keyserver hkp://keyserver.ubuntu.com:80 422C4D99 \
    && apt-get update \
    && apt-get install -y chaste-dependencies \
    && apt-get clean \  
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# See https://github.com/phusion/baseimage-docker/issues/186
RUN touch /etc/service/syslog-forwarder/down

# The entrypoint script below will ensure our new chaste user (for doing builds)
# has the same userid as the host user owning the source code volume, to avoid
# permission issues.
# Based on https://denibertovic.com/posts/handling-permissions-with-docker-volumes/
COPY entrypoint.sh /usr/local/bin/entrypoint.sh

# Hook to link to host chaste source folder, and set it as the working dir
VOLUME /usr/src/chaste
WORKDIR /usr/src/chaste

# Use baseimage-docker's init system, and switch to the chaste user running
# bash as a login shell by default (see entrypoint.sh).
# If no specific command is given the default CMD will drop us into an
# interactive shell.
ENTRYPOINT ["/sbin/my_init", "--quiet", "--", "/usr/local/bin/entrypoint.sh"]
CMD ["bash -i"]
