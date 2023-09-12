# syntax=docker/dockerfile:1

ARG BASE_TAG=latest
FROM chaste/base:${BASE_TAG}
# FROM chaste/develop

# Build Chaste: TAG can be a branch or release ('-' skips by default)
ARG GIT_TAG=.
ENV GIT_TAG=${GIT_TAG}
ADD . ${CHASTE_SOURCE_DIR}
RUN git config --global --add safe.directory /home/chaste/src
RUN echo "Building Chaste ${GIT_TAG} in $PWD..." && ls
RUN build_chaste.sh ${GIT_TAG}

# Automatically mount the home directory in a volume to persist changes made there.
# NOTE: After declaring the volume, changes to the contents during build will not persist.
VOLUME "${CHASTE_DIR}"

# Optionally run a test suite before finalising the image.
# NOTE: These test outputs will not appear in the volume. 
ARG TEST_SUITE=-
ENV TEST_SUITE=${TEST_SUITE}
RUN test.sh ${TEST_SUITE}
