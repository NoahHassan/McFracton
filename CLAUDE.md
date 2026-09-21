# McFracton — working notes for Claude

Lattice Monte Carlo code (C++20). `McMachine` runs Metropolis on an abstract `System`; the
`gui/` app visualises a running simulation with SFML 2.6 + Dear ImGui.

`REFACTOR_PLAN.md` is the agreed architecture cleanup plan. Read it before making structural
changes, and keep it up to date when decisions change.

## Rules
1. **Never `git commit` without asking first.**
2. **Do not change the physics.** Energies, update proposals, defect definitions and the order of
   floating-point operations stay as they are. Physics problems that are noticed get written into
   Appendix A of `REFACTOR_PLAN.md` for Noah, never silently fixed.
3. Visual Studio stays a first-class way to work: `McFracton.sln` must keep building and running.
   Every file added or removed has to be reflected in `McFracton.vcxproj` and `.vcxproj.filters`.
4. Changes that alter logged numbers need agreement first, and belong before the regression
   baseline is captured.

## Layout
```
McFracton/core/include, core/src   physics + Monte Carlo core, no SFML
McFracton/gui                      SFML/ImGui app
McFracton/tests/golden             regression baselines
```

## Build
Visual Studio (primary, Windows):
```sh
"/c/Program Files/Microsoft Visual Studio/18/Community/MSBuild/Current/Bin/MSBuild.exe" \
    McFracton.sln -p:Configuration=Release -p:Platform=x64 -v:minimal -nologo
```
Output: `x64/Release/McFracton.exe`. Dependencies come from vcpkg at `C:/src/vcpkg`
(sfml 2.6.1, imgui, imgui-sfml), via the user-wide vcpkg MSBuild integration.

CMake (for the Linux/Slurm cluster; also usable on Windows):
```sh
CMAKE="/c/Program Files/Microsoft Visual Studio/18/Community/Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/bin/cmake.exe"
"$CMAKE" --preset windows-release && "$CMAKE" --build --preset windows-release
```
Output: `build/`. Sources are globbed, so new files need no CMake edit.

## Verification
Both builds must pass at every checkpoint. The regression test compares against the golden files:
```sh
build/windows-release/mcf_run --regression --compare McFracton/tests/golden
```
It must be bit-identical unless a change to logged numbers was explicitly agreed.

## Notes
- The GUI runs on SFML 2.x API (`sf::Event` polling, `sf::VideoMode(w, h)`), not SFML 3.
- `QXYSquare` is out of scope: it must keep compiling, but it is excluded from the regression test.
- `x64/` and `build/` are build output and are gitignored.
