# McFracton — Architecture Cleanup Plan

## Context

McFracton is a C++ lattice Monte Carlo code: `McMachine` runs Metropolis on an abstract `System`, with
5 implementations (classical XY, quantum XY, U(1) gauge 2+1D / 3+1D, Spiderweb rank-2 fracton) plus an
SFML/imgui GUI. As the project grew, the abstractions stopped fitting:

- `System` carries a fixed `Observables` struct (helicity, Polyakov loop, defects a/b …) most systems don't use;
  `Measure(double T)` forces temperature on everyone; the `site`/`plaq` split is vestigial (the plaq side is never
  used by `McMachine`, yet every system implements `proposePlaqFlip`/`UpdatePlaq`/`plaq_fields`).
- Every system re-implements lattice indexing, periodic wrap, `mapToCircle`, `sgn`, `#define PI`.
- `McMachine` has a hard-coded log header, and opens its log file even when only the GUI runs.
- The GUI has three near-identical canvases (`Canvas`, `HyperCanvas`, `SpiderCanvas`), each hard-typed to one
  system; `Main.cpp` must be edited to switch systems.
- Only a Visual Studio project exists, but the core must run headless on a **Slurm cluster (Linux)**.

**Goal:** a system-agnostic architecture — generic observables, one generic GUI canvas, a headless CLI,
CMake build — **without changing any physics**, and **without blowing up the file count**.

### Ground rules (agreed)
1. Never `git commit` without asking first.
2. Do not change physics. Formulas keep their exact order of operations, verified by bit-identical regression.
   Physics bugs found along the way are **only listed** (Appendix A), not fixed.

### Decisions made
- Depth: change the interface **and** dedupe helpers inside the systems.
- Build: **Visual Studio stays the primary build on Windows.** `McFracton.sln` / `McFracton.vcxproj` are kept and
  keep working exactly as now — open, edit, build, run and play the GUI from VS, with output in `x64/` as before.
  A `CMakeLists.txt` is added **alongside** it for the Linux cluster, building the core and the headless CLI
  (the GUI is optional and off by default there). CMake builds into `build/`, so the two never collide.
  To stop the two file lists from drifting, CMake picks up sources with `file(GLOB CONFIGURE_DEPENDS)`, so adding
  a file in Visual Studio is automatically picked up by the cluster build with no CMake edit.
  A second small project, `McFractonRun.vcxproj` (the headless CLI), is added to the same solution, so the
  regression test and batch runs are also runnable from inside VS. I keep both projects up to date whenever I add
  a file; if you ever only want the CLI on the cluster, that second project can simply be dropped.
- Observables become generic doubles, so the defect-count mean/variance lose their integer truncation
  (an intended change, in those columns only).
- `mapToCircle`: only the **period-2** variant (the one now in every system) goes into `MathUtil.h`.
- **No `Annealer`, `Statistics` or `ResultsWriter` files.** The temperature schedule, the mean/variance/
  autocorrelation helpers and the TSV logging all stay inside `McMachine`, as they are today.
- Overrelaxation keeps the base-class **throw** as its guard, so an unimplemented overrelaxation fails loudly.
  The GUI loses its overrelax checkbox entirely.
- `measure(double temperature)` keeps its explicit temperature argument, as in the original. Only the fixed
  `Observables` struct goes away; what changes is that each system declares *which* observables it reports.
- **QXYSquare is out of scope**: it stays in the tree, compiles, but is excluded from the regression tests.

---

## Target architecture

```
McFracton/                      (repo root)
├─ McFracton.sln                KEPT — primary Windows build (GUI + CLI projects)
├─ CMakeLists.txt               added for the cluster; option(MCF_BUILD_GUI, OFF on Linux)
├─ CMakePresets.json            linux-release (no GUI), windows-release (vcpkg toolchain, optional)
├─ CLAUDE.md                    rules + build/test commands
├─ REFACTOR_PLAN.md             copy of this plan
├─ scripts/slurm_array.sh       example array job
└─ McFracton/
   ├─ McFracton.vcxproj         KEPT — the GUI app, as today
   ├─ McFractonRun.vcxproj      NEW — the headless CLI, same solution
   ├─ core/include/  core/src/          → mcf_core (no SFML)
   │   System.h         abstract system (new interface below)
   │   McMachine.h      update kernel + annealing schedule + statistics + TSV logging (as today)
   │   Lattice.h        NEW: PeriodicLattice — extents, index(coords), coords(index), neighbor(i, axis, ±k)
   │   MathUtil.h       NEW: kPi (= std::numbers::pi, the same double as the old macro), sgn, mapToCircle
   │   SystemRegistry.h NEW: name → {parameter spec, factory}; the one thing that lets a new system reach
   │                    both the CLI and the GUI without either of them being edited
   │   systems/         XYSquare, QXYSquare, AbelianGaugeSquare, AbelianGaugeCube, Spiderweb
   ├─ cli/main.cpp      → mcf_run   (headless run; also `--regression` mode; key=value config parsing lives here)
   └─ gui/                          → mcf_gui (SFML 2.6 + imgui-sfml via vcpkg)
       Main.cpp         system picker (registry combo + params + Rebuild), controls, energy plot
       LatticeView.h/.cpp  the ONE generic canvas
       ColorMaps.h      GreenRedUniform, NormalMapYellow, BlackWhite, RedWhiteBlue (moved verbatim)
       BufferedArray.h  (moved from core; GUI-only)
```
Net new files: `Lattice.h`, `MathUtil.h`, `SystemRegistry.h` in core; in the GUI, three canvases collapse into one.

### New `System` interface (sketch)
```cpp
struct LatticeInfo { std::vector<int> extents; std::vector<std::string> axis_names; }; // {L,L,Nt}, {"x","y","t"}
enum class ChannelKind { Angle, Signed, Integer, Magnitude };   // picks a GUI colormap
struct Channel { std::string name; ChannelKind kind; };

class System {
public:
    explicit System(int n_variables);
    virtual ~System() = default;
    int    numVariables() const;
    double variable(int i) const;

    // Metropolis
    virtual double getEnergy() const = 0;
    virtual double proposeUpdate(int i, double delta) const = 0;   // was proposeSiteFlip
    void           applyUpdate(int i, double delta);                 // was UpdateSite (identical everywhere)
    virtual void   overrelax(int i) { throw std::logic_error("overrelax not implemented"); }

    // Measurement: names are declared once, measure() returns values in the same order
    virtual std::vector<std::string> observableNames() const = 0;
    virtual std::vector<double>      measure(double temperature) const = 0;

    // Visualisation as pure data, no SFML: named scalar fields on the lattice
    virtual LatticeInfo          lattice()  const = 0;
    virtual std::vector<Channel> channels() const = 0;
    virtual void fillChannel(int channel, std::vector<double>& out) const = 0;  // one value per lattice site
protected:
    std::vector<double> fields;
};
```
The overrelax guard stays exactly as it is in spirit; only the thrown type changes from a `const char*` to
`std::logic_error`, so it can be caught and printed with `.what()`. Say the word and it stays a bare `throw("…")`.

Removed: `Observables` struct, `plaq_fields`, `n_plaq_variables`, `proposePlaqFlip`, `UpdatePlaq`, `getPlaq`,
`getLocalEnergies` (it becomes a Spiderweb channel), `QXYSquare::LogToFile`.
Observable names keep the old log-header spellings (`Energy`, `Helicity Modulus`, `defects_a`, `Polyakov Loop`),
so analysis scripts can find columns by name; each system now logs only the columns it actually computes.

### Generic GUI canvas (`LatticeView`)
- Reads `lattice()`: you pick which two axes are shown (default the first two), and every remaining axis gets a
  slice slider. One view covers 2D, 2+1D and 3+1D, replacing `Canvas`, `HyperCanvas` and `SpiderCanvas`.
- The channel combo comes from `channels()`; the colormap follows `ChannelKind` (Angle→NormalMapYellow,
  Integer/Signed→RedWhiteBlue, Magnitude→BlackWhite with autoscale). An optional overlay channel draws only
  nonzero cells, replacing today's "draw monopoles on top".
- Cells are drawn as one `sf::VertexArray` (quads with a 1px gap for the grid look) instead of a
  `shared_ptr<sf::Shape>` per cell, which makes `Shape`, `Square` and `Vec2D` unnecessary.
- Keys: ←/→ step the first slice axis, Enter steps the second, Space pauses — as today.
- The overrelax checkbox **and** the `if (params.overrelax) machine.Overrelax(500);` call are removed from the
  render loop. `NumericalParams::overrelax` stays, because the headless annealing run still uses it.
- Adding a system then needs **zero GUI changes**: register it and implement `lattice/channels/fillChannel`.

| System | Channels (reusing existing functions unchanged) |
|---|---|
| XYSquare | θ (Angle); vortices (Integer, from `getVortices`) |
| QXYSquare | θ (Angle) |
| AbelianGaugeSquare | A_x, A_y, A_t (Angle); flux_xy (`getFluxes_z`); monopoles (`getMonopoles`) |
| AbelianGaugeCube | A_x, A_y, A_z, A_t; flux; monopoles |
| Spiderweb | A_0, A_xx, A_xy (Angle); local energy (`getLocalEnergies`) |

---

## Phases
Each phase ends with a build, the regression check, and a checkpoint where I ask whether to commit.

### Phase 0: Environment and safety net (no refactor yet)
1. **Repair the build.** The working tree does not compile or link right now, after the physics edits:
   - `AbelianGaugeCube::Measure` still calls the deleted `getMonopoles()`; remove that defect-count block
     (`n_defects_a/b` stay 0 for the Cube) and drop the declaration from `AbelianGaugeCube.h`.
   - `gui/HyperCanvas.h` calls `field.getMonopoles()`; remove its monopole path (that file disappears in Phase 5).
   - `Spiderweb.h` still declares `OverrelaxSite` as an override while the definition is commented out; remove the
     declaration, so Spiderweb falls back to the base class throw — which is exactly the intended guard.
2. Copy this plan to `REFACTOR_PLAN.md`, write `CLAUDE.md` with the rules and build commands, save the rules to memory.
3. Build setup, without disturbing your Visual Studio workflow:
   - `McFracton.vcxproj` keeps building the GUI exactly as now; I add every new file to it (and to `.filters`)
     as I go, so the solution never breaks.
   - `CMakeLists.txt` + `CMakePresets.json` build the **current** code with globbed sources, verified locally with
     the VS-bundled CMake/Ninja and the vcpkg toolchain at `C:/src/vcpkg`.
   - `.gitignore`: add `build/`, untrack `imgui.ini`, ignore run outputs (`results/`, `*.tsv`). `x64/` is already
     ignored and stays where it is.
4. Minimal **seed plumbing**: an optional `seed` parameter for `McMachine` and each system's RNG, defaulting to
   `random_device`, so behavior is unchanged.
5. **Deliberate numerical changes, all made before the baseline is captured** (each one changes logged numbers,
   which is why they belong here and not in a later phase):
   - Drop the `(float)` cast in `Measure`'s call to `Sweep`, and make the acceptance RNG `double` rather than
     `float`. Both are pure precision losses with nothing gained.
   - Fix `XYSquare::getPlaqConnectedCluster`: it currently uses open-boundary `(size-1)` indexing while being
     called for `size*size` plaquettes, which reads out of range. It becomes periodic `size` indexing, matching
     the rest of the system.
   - `Timer::elapsed()` returns seconds (`std::chrono::duration<float>`) instead of raw ticks. The CLI then uses
     it to report wall time per temperature.
   - **Freeze δ during measurement**: it adapts freely while thermalizing, then stays fixed for that
     temperature's measurement sweeps, which removes the adaptive-chain technicality (Appendix A.3).
   - **Error columns become the standard error** `s/√N`, with `s²` using the `N−1` convention, instead of the
     raw variance (Appendix A.4). Column count and header names stay as they are.
6. Portability fixes that cannot change results: missing includes (`<cmath>` in McMachine.cpp, `<vector>` in
   System.h), and the default argument moved from McMachine.cpp into the header.
7. `mcf_run --regression`: for each system except QXYSquare (small L, fixed seed) run N sweeps and a few `Measure`
   calls at 2–3 temperatures, dumping energy and all observables at `%.17g`; golden files in
   `McFracton/tests/golden/`. **This is captured last**, once steps 1–6 are done and you confirm the physics is
   settled, so every later "bit-identical" claim is anchored to a state you are no longer editing.

### Phase 1: Shared utilities and dedupe inside the systems
- Add `MathUtil.h` and `Lattice.h`; replace each system's `to_site_index`/`index_from_site`/modular wraps/
  `mapToCircle`/`sgn`/`PI`. Since all five now share one `mapToCircle`, this is a pure move.
- Spiderweb: `getEnergy()` and `getEnergy(vector&)` both call one `computeLocalEnergies()` helper with the same
  summation order.
- Remove dead code: commented-out blocks, `ColorVortices`/`meh`, unused `direction` variables, the unused
  `XYSquare(size, temperature)` parameter, `const` on by-value returns.
- Regression: **bit-identical**.

### Phase 2: New `System` interface and generic observables
- Apply the interface above; drop the plaquette machinery; move each old `Measure(T)` body into
  `observableNames`/`measure`, with XYSquare's helicity modulus as its one derived observable.
- `McMachine` builds its log header from the observable names, and only opens the log file for an actual
  annealing run, so the GUI stops creating an empty log file at every launch. Schedule, statistics and writing
  stay inside `McMachine`.
- `QXYSquare::measure` returns `{Energy}` instead of throwing — an intentional improvement, which makes a
  headless run of that system possible for the first time.
- Today's numerical quirks stay exactly as they are (`(float)temperature` cast, float acceptance RNG, delta
  adaptation during measurement, the `max_therm_sweeps` clamp). They are listed in Appendix A.
- Log format: `#`-prefixed provenance lines (git hash, parameters, seed) come first, and the column-header line
  stays the first non-comment line, so parsers reading with `comment='#'` keep working.
- Regression: bit-identical, except the defect mean/variance columns, which are now doubles, as agreed.

### Phase 3: RNG ownership (the one step with a re-baseline)
- Systems stop owning RNGs. `McMachine` owns one seeded `std::mt19937` and passes it to `system.randomize(rng)`
  and `system.overrelax(i, rng)`, which makes a run reproducible from a single seed.
- That changes the random stream, so verification here is **statistical**: a longer seeded run at a few T per
  system, with mean energy and observables agreeing within error bars. Then new golden files.
- Note: the same seed gives different streams on MSVC and GCC, because `std::uniform_*_distribution` is
  implementation-defined. Runs are reproducible per platform.

### Phase 4: Headless CLI for Slurm
- `SystemRegistry` plus `mcf_run`: a key=value config file with `--key=value` overrides for the system, sizes,
  couplings, `NumericalParams`, `seed` and `out`. Config parsing lives in `cli/main.cpp`, with no new dependency.
- `scripts/slurm_array.sh` as an example, and a Linux build check if you can run
  `cmake --preset linux-release && cmake --build --preset linux-release` on the cluster (I can't reach it).

### Phase 5: Generic GUI
- `LatticeView`, `ColorMaps.h`, registry-driven `Main.cpp`. Delete `Canvas.h`, `HyperCanvas.h`, `SpiderCanvas.h`,
  `Shape.*`, `Square.h`, `Vec2D.h`, `Timer.h` (unused).
- Manual check: every system renders every channel, slicing works along all axes, and the energy plot and pause work.

### Phase 6: Documentation
- `README.md`: how to build and run both ways — open `McFracton.sln` in Visual Studio for the GUI, and
  `cmake --preset linux-release` on the cluster for the headless runs. Nothing is deleted; the solution, the
  project files and `x64/` all stay.

---

## Verification
- **Build, both ways, at every checkpoint**: MSBuild on `McFracton.sln` (so your Visual Studio workflow is proven
  to still work, output in `x64/`), and `cmake --build --preset windows-release` with the VS-bundled `cmake.exe`
  (output in `build/`), which stands in for the cluster build.
- **Regression**: `mcf_run --regression --compare McFracton/tests/golden`. Bit-identical through Phases 0–2
  (except the agreed defect columns); statistical agreement after Phase 3, then new golden files.
- **CLI smoke test**: a short `mcf_run` per system produces a well-formed TSV.
- **GUI**: launch `mcf_gui` and check by hand.
- Checkpoint after each phase: I summarize the diff and ask before any commit.

---

## Appendix A: Physics and numerics issues

### Resolved by you (2026-09-21)
`mapToCircle` unified to the period-2 variant everywhere; `cos(2πx)` → `cos(πx)`; XYSquare's vortex threshold
raised to ±2; the Cube's `getMonopoles` deleted as unphysical and its Polyakov loop fixed to direction 3 with the
`PI` factor; `getFluxes_z` fixed to `site*6 + type`; `Spiderweb::get_field` fixed to `site*3 + type` and its
`OverrelaxSite` commented out; `current_measurement_sweeps` now clamped by `max_measure_sweeps`.
The `+`/`−` signs in the Spiderweb Hamiltonian are intended, and gauge "overrelaxation" keeps its name.

### Still open, for you (I will not touch these)
1. **`XYSquare::getSinSqrX` (line 60) still uses `sin(2.0 * PI · Δθ)`** while the energy is now `cos(PI · Δθ)`.
   The helicity modulus is the twist derivative of the energy, so its `sin` should carry the same argument as the
   `cos` it comes from, and the `π` prefactors in the formula change with the convention. This is the last
   `2.0 * PI` left in the core.
2. QXYSquare is out of scope, and its `Measure` still throws. It is excluded from the regression tests and is left
   untouched; say the word if you'd rather delete it.

### Numerics being changed in Phase 0 (agreed)
3. **Adaptive step size δ.** Tuning δ toward ~50% acceptance is right, and the update rule has its fixed point
   exactly at 0.5. The subtlety is only *when* it is tuned: a Metropolis proof assumes a fixed proposal kernel, so
   adapting δ from the chain's own history during the measurement phase makes the chain adaptive rather than a
   homogeneous Markov chain. → **δ is frozen for the measurement sweeps.**
4. **The `D*` columns are variances of single measurements**, `s² = (1/N)Σ(x−x̄)²`, not uncertainties of the
   mean. As an error bar the mean needs `s/√N` (with `s²` using `N−1`). Since measurements are taken every
   ~20·τ_int sweeps, they are effectively independent, so no extra τ factor is needed.
   → **The columns become `s/√N`**, keeping their names and positions.
