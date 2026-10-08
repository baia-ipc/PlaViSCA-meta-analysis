#!/usr/bin/env bash
# Idempotently ensures Bioconductor "data-only" packages are actually
# installed into the pixi environment's R library.
#
# Root cause: bioconda's Bioconductor *data* recipes (e.g.
# bioconductor-genomeinfodbdata) ship the R package metadata via conda but
# fetch the actual (large) data tarball at install time through a conda
# `post-link.sh` hook (see the packaged
# .pixi/envs/default/bin/.bioconductor-genomeinfodbdata-post-link.sh and
# installBiocDataPackage.sh). pixi/rattler deliberately do not execute
# conda post-link/pre-unlink scripts (unlike conda/mamba), so the data
# payload never gets downloaded - `library(GenomeInfoDbData)` (and
# therefore GenomeInfoDb/rtracklayer/DropletUtils, which depend on it)
# fails even though conda-meta reports the package as installed.
#
# This script reproduces exactly what that post-link.sh would have done,
# but as an explicit, version-controlled, idempotent step run from the
# pixi environment's own activation, targeting only the controlled pixi
# env's R library (never a personal/user library).
#
# IMPORTANT: pixi activation scripts are SOURCED into the activating shell,
# not executed as an independent subprocess - this function therefore never
# calls `exit`/`set -e` at the top level (either would terminate/corrupt the
# calling activation shell and silently abort the rest of PATH/env setup,
# falling back to the system R). All logic is scoped inside a function.

_plavisca_ensure_genomeinfodbdata() {
  local prefix="${CONDA_PREFIX:-}"
  if [ -z "$prefix" ] || [ ! -x "$prefix/bin/R" ]; then
    return 0
  fi

  local rlib="$prefix/lib/R/library"
  if [ -d "$rlib/GenomeInfoDbData" ]; then
    return 0
  fi

  echo "[ensure_bioc_data_packages] GenomeInfoDbData missing from $rlib (pixi does not run conda post-link.sh) - installing..." >&2

  local pkg_version="1.2.11"
  local pkg_md5="2a4cbfc2031992fed3c9445f450890a2"
  local pkg_fn="GenomeInfoDbData_${pkg_version}.tar.gz"
  local urls=(
    "https://bioconductor.org/packages/3.18/data/annotation/src/contrib/${pkg_fn}"
    "https://bioarchive.galaxyproject.org/${pkg_fn}"
    "https://depot.galaxyproject.org/software/bioconductor-genomeinfodbdata/bioconductor-genomeinfodbdata_${pkg_version}_src_all.tar.gz"
  )

  local staging="$prefix/share/genomeinfodbdata-${pkg_version}"
  mkdir -p "$staging"
  local tarball="$staging/$pkg_fn"

  local success=0
  local url
  for url in "${urls[@]}"; do
    if curl -fsSL --max-time 120 "$url" -o "$tarball" && echo "${pkg_md5}  ${tarball}" | md5sum -c - >/dev/null 2>&1; then
      success=1
      break
    fi
  done

  if [ "$success" != 1 ]; then
    echo "[ensure_bioc_data_packages] ERROR: could not download/verify $pkg_fn from any mirror; library(GenomeInfoDbData) will fail" >&2
    rm -rf "$staging"
    return 1
  fi

  "$prefix/bin/R" CMD INSTALL --library="$rlib" "$tarball" >&2
  rm -rf "$staging"
  echo "[ensure_bioc_data_packages] GenomeInfoDbData installed successfully." >&2
  return 0
}

_plavisca_ensure_genomeinfodbdata || echo "[ensure_bioc_data_packages] continuing despite install failure - some R packages may not load" >&2
unset -f _plavisca_ensure_genomeinfodbdata
