# Interactive Earth

<img src="cover.png" align="center" width="512">

Interactive Earth is a teaching tool for demonstrating mantle convection and its dependence on Rayleigh number.
It enables a mantle convection model to be run much more quickly than the codes typically used in research.
It solves the same governing equations, but makes comprimises for accuracy for the sake of achieving speed.
This speed enables the code to be very responsive to user input.

The simulations can be run in a browser without running any local code [here](http://ian-r-rose.github.io/interactive_earth/thermal_hires.html).

More about the project can be found [here](http://ian-r-rose.github.io/interactive_earth/about.html).

An introduction to the underlying science can be found [here](http://ian-r-rose.github.io/interactive_earth/explanation.html).

## Building

The build uses CMake. The easiest way to get a known-good toolchain
(compiler, CMake, Ninja, FFTW, and the X11/OpenGL libraries) is [pixi](https://pixi.sh),
which installs pinned versions from `pixi.lock` into the project directory:

```
pixi run start    # build and run the interactive app
pixi run test     # physics regression tests
pixi run bench    # headless solver benchmark, e.g. pixi run bench --ra 1e6 --threads 1
```

Windowing and input use the vendored [sokol_app.h](https://github.com/floooh/sokol) in `third_party/`.
Without pixi, install those dependencies however you like and use the presets directly:

```
cmake --preset native-release
cmake --build --preset native-release
ctest --preset native-release
./build/native-release/interactive_earth
```

There is also a `native-debug` preset with AddressSanitizer and UndefinedBehaviorSanitizer.

## Building for the web

The WebAssembly build uses the Emscripten SDK version pinned in `.emsdk-version`,
which `pixi run setup-emsdk` installs into `.emsdk/`. FFTW is downloaded and built
from source with the same toolchain.

```
pixi run setup-emsdk   # once
pixi run build-wasm
pixi run serve-wasm    # then open http://127.0.0.1:8000
pixi run test-wasm     # run the physics tests under node
pixi run bench-wasm    # single-threaded benchmark under node
```

`web/dev.html` is a minimal page for local testing; add `?log` to the URL to also
log the timing overlay to the console, and `?grid=` to set the grid size (see below),
e.g. `http://127.0.0.1:8000/?grid=2048x256`.

## Resolution and run-time settings

The grid size sets both the cost of a step and the highest Rayleigh number the model
allows, which goes as the cube of the number of radial cells (so that the thermal
boundary layers stay resolved). The default is 1536x192, for a maximum Ra of about 1.4e8.
To choose another, set `IE_GRID` (or `?grid=` on the web page):

```
IE_GRID=2048x256 pixi run start   # maximum Ra about 3.4e8
IE_GRID=1024x128 pixi run start   # maximum Ra about 4.2e7, and faster
```

The app takes as many solver steps per frame as fit in a 60 fps frame, and reports the
step cost, steps per frame, and frame rate in the terminal (or the web page's overlay)
about once a second. Natively it uses OpenMP threads; set `OMP_NUM_THREADS` to limit them.

For regression checks, `IE_SCREENSHOT=out.ppm` seeds convection with fixed heat sources,
runs `IE_SCREENSHOT_STEPS` steps (default 200), saves the window to that file, and quits.
User input is ignored in this mode, so the result is repeatable.

## Development

Code is formatted with clang-format (`.clang-format`, pinned to one major version by pixi):

```
pixi run format          # format all sources
pixi run install-hooks   # once per clone: format staged files on every commit
pixi run lint            # clang-tidy (.clang-tidy)
```

CI (`.github/workflows/ci.yml`) checks formatting, lints, builds, and runs the tests for
both the native and WebAssembly builds, on pushes and pull requests to `main`.

## User input
With the simulation running, the user can make the following changes:
- put in positive temperature anomalies (add heat) by clicking
- put in negative thermal anomalies (remove heat) by right-clicking (or Ctrl-clicking)
- change the Rayleigh number (the vigor of convection) by scrolling
- press space to pause convection and switch to seismic mode, where clicking makes
  earthquakes and right-clicking places the seismometer; press space again to resume
- press `c` to cycle through the color scales
- press Backspace to turn advection and diffusion off (and on again), to paint with heat
- press Tab to switch between the temperature and composition fields, when the model
  includes a chemical field; clicking then adds composition
- press Escape to quit (natively)

## Changing simulation behavior for advanced users
The behavior of the simulation can be modified by changing parameters within the upper part of `src/main.cpp`. The parameters that readily exposed for being modified include:
- the default resolution of the simulation (and thereby how quickly it runs), which
  `IE_GRID` overrides
- whether or not to do a true polar wander (TPW) calculation
- whether to included a chemical field
- whether to render the simulation with flattened with an equatorial bulge
