#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# System libraries for Code Ocean / Ubuntu 22.04
# ============================================================
# This script installs system libraries only.
# It intentionally does not install Ubuntu r-cran-* or r-bioc-* packages,
# because mixing Ubuntu R binary packages with the CRAN R runtime used by
# Code Ocean can cause version conflicts.
#
# Do not run add-apt-repository in the capsule unless required.
# ============================================================

apt-get update

apt-get install -y \
  build-essential \
  gfortran \
  pkg-config \
  cmake \
  git \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev \
  libhdf5-dev \
  libglpk-dev \
  libmagick++-dev \
  jags \
  libgdal-dev \
  libgeos-dev \
  libproj-dev \
  libudunits2-dev

echo "System dependency installation finished."
