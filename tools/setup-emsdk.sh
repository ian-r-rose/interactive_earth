#!/usr/bin/env bash
# Install the Emscripten SDK version pinned in .emsdk-version into .emsdk/
# (git-ignored). emsdk's tags match SDK versions, so the installer script
# itself is pinned too.
set -euo pipefail
cd "$(dirname "$0")/.."
version=$(cat .emsdk-version)
if [ ! -d .emsdk ]; then
  git clone --quiet https://github.com/emscripten-core/emsdk.git .emsdk
fi
git -C .emsdk fetch --quiet --tags
git -C .emsdk -c advice.detachedHead=false checkout --quiet "$version"
.emsdk/emsdk install "$version"
.emsdk/emsdk activate "$version"
