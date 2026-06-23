#!/bin/bash
# Print the conda package version for the Chaste source tree at $1 (default: .).
src="${1:-.}"
tag=$(git -C "${src}" tag --points-at HEAD 2>/dev/null | head -1)
if [ -n "${tag}" ]; then
  echo "${tag#v}"
else
  short_hash=$(git -C "${src}" rev-parse --short=8 HEAD)
  echo "0.dev.g${short_hash}"
fi
