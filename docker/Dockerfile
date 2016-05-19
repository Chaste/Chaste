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

#COPY set-home-permissions.sh /etc/my_init.d/set-home-permissions.sh
#RUN chmod +x /etc/my_init.d/set-home-permissions.sh

# Set up chaste user so that we do not build as root
RUN useradd -m -s /bin/bash -G sudo,docker_env chaste && \
    echo "chaste:docker" | chpasswd && \
    echo "chaste ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers
USER chaste

# Hook to link to host chaste source folder, and set it as the working dir
VOLUME /usr/src/chaste
WORKDIR /usr/src/chaste

# Use baseimage-docker's init system, and switch to the chaste user running
# bash by default.  The -l makes bash behave like a login shell, so the USER
# environment variable will be set.  The -c at the end of the ENTRYPOINT
# means we can do things like "docker run -it ... scons" to build Chaste.  If
# no specific command is given the default CMD will drop us into an
# interactive shell.
ENTRYPOINT ["sudo", "/sbin/my_init", "--quiet", "--", "sudo", "-u", "chaste", "/bin/bash", "-l", "-c"]
CMD ["bash -i"]
