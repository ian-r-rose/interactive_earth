#!/usr/bin/env bash
# Run a command with the pinned Emscripten SDK on the PATH, e.g.
#   tools/with-emsdk.sh emcmake cmake --preset wasm-release
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
if [ ! -f "$root/.emsdk/emsdk_env.sh" ]; then
  echo "Emscripten SDK not installed; run tools/setup-emsdk.sh (or pixi run setup-emsdk)" >&2
  exit 1
fi
source "$root/.emsdk/emsdk_env.sh" > /dev/null 2>&1
exec "$@"
