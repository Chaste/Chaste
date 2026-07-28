#!/bin/bash
# Print the conda package version for the Chaste source tree at $1 (default: .).
src="${1:-.}"
tag=$(git -C "${src}" tag --points-at HEAD 2>/dev/null | head -1)
if [ -n "${tag}" ]; then
  echo "${tag#v}"
else
  latest_tag=$(git -C "${src}" describe --tags --abbrev=0 2>/dev/null)
  latest_tag="${latest_tag:-0}"
  count=$(git -C "${src}" rev-list --count "${latest_tag}..HEAD" 2>/dev/null || git -C "${src}" rev-list --count HEAD)
  short_hash=$(git -C "${src}" rev-parse --short=8 HEAD)
  echo "${latest_tag#v}.dev${count}.g${short_hash}"  # PEP 440 dev version: LATEST_TAG.devN.gHASH
fi
