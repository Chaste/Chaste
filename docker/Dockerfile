# syntax=docker/dockerfile:1

ARG BASE_TAG=latest
FROM chaste/base:${BASE_TAG}

# Build Chaste: GIT_TAG can be: a branch or release; '-' to skip; '.' to build local code
ADD . ${CHASTE_SOURCE_DIR}
# Fix git permissions issue CVE-2022-24765
RUN git config --global --add safe.directory "${CHASTE_SOURCE_DIR}"
ARG GIT_TAG=.
ENV GIT_TAG=${GIT_TAG}
RUN build_chaste.sh ${GIT_TAG}

# Automatically mount the home directory in a volume to persist changes made there.
# NOTE: After declaring the volume, changes to the contents during build will not persist.
VOLUME "${CHASTE_DIR}"

# Optionally run a test suite before finalising the image.
# NOTE: These test outputs will not appear in the volume. 
ARG TEST_SUITE=-
ENV TEST_SUITE=${TEST_SUITE}
RUN test.sh ${TEST_SUITE}
