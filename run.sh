#!/usr/bin/env bash
# Run pipeline stages in the pixi environment and keep a log of each.
#
# Usage:  ./run.sh all            run every stage, in order
#         ./run.sh 30 40          run the given stages
#         ./run.sh list           list the stages
#         ./run.sh check-inputs   verify the input files against their checksums
#
# Logs go to <workspace>/runs/<run>/logs/<stage>.log (see config/local.env).
set -euo pipefail
cd "$(dirname "$0")"

# settings from config/local.env, unless already set in the environment
if [ -f config/local.env ]; then
  while IFS='=' read -r key value; do
    [[ "$key" =~ ^[A-Z_]+$ ]] && [ -z "${!key:-}" ] && export "$key=$value"
  done < <(grep -E '^[A-Z_]+=' config/local.env)
fi
: "${PLAVISCA_WORKSPACE:?set it in config/local.env (see config/local.env.example)}"
log_dir="$PLAVISCA_WORKSPACE/runs/${PLAVISCA_RUN:-main}/logs"
mkdir -p "$log_dir"

stages=$(ls R/[0-9][0-9]_*.R | xargs -n1 basename | sed 's/\.R$//')

if [ "${1:-}" = "check-inputs" ]; then
  inputs="$PLAVISCA_WORKSPACE/inputs"
  (cd "$inputs" && tail -n +2 "$OLDPWD/metadata/input_files.tsv" | awk -F'\t' '$3 ~ /^[0-9a-f]{64}$/ {print $3 "  " $1}' | sha256sum -c --quiet) &&
  (cd "$inputs/starsolo_counts" && sha256sum -c --quiet "$OLDPWD/metadata/starsolo_counts.sha256") &&
  (cd "$inputs/study_count_objects" && sha256sum -c --quiet "$OLDPWD/metadata/study_count_objects.sha256") &&
  echo "all input files match their checksums"
  exit
fi

if [ "${1:-}" = "" ] || [ "$1" = "list" ]; then
  echo "$stages"; exit 0
fi
[ "$1" = "all" ] && set -- $(echo "$stages" | cut -c1-2)

for n in "$@"; do
  stage=$(echo "$stages" | grep "^$n" || true)
  [ -n "$stage" ] || { echo "no stage $n" >&2; exit 1; }
  echo "=== $stage ($(date '+%F %T'))"
  pixi run Rscript --vanilla "R/$stage.R" 2>&1 | tee "$log_dir/$stage.log"
done
