#!/usr/bin/env bash
# Backward-compatible entrypoint. Use start_nagini.sh for new runs.
set -euo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "$script_dir/start_nagini.sh" "$@"
