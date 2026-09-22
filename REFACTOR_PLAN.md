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
├─ config/spiderweb.cfg         example run configuration
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

**As built (deviations from the sketch above, all agreed after Phase 1 landed):**
- `LatticeInfo` and the virtual `lattice()` are **not** used. Phase 1's `mcf::PeriodicLattice` already
  carries extents and axis names, so a second struct would have duplicated it: `System` now *holds* the
  lattice as a protected member, hands it out through a non-virtual `getLattice()`, and its constructor
  takes `(PeriodicLattice, n_variables)` and sizes `fields` itself. That also removed the `lattice`
  member and the `fields` allocation from all five systems.
- Error columns are named `D` + the observable name (`DEnergy`, `DHelicity Modulus`), since the old
  hand-written abbreviations (`DE`, `DHM`, `Dna`) cannot be derived from a name. The value columns keep
  the old spellings exactly, which is what analysis scripts key on.
- The `#` provenance lines carry the resolved seed (so a `seed = 0` run records the seed it actually
  drew), the `NumericalParams` and the variable count. **The git hash is deferred to Phase 4**, where
  the CLI can pass it in; embedding it now would have needed a generated header in both build systems.
- `getLocalEnergies()` left the base interface but stays a `Spiderweb` method, as the source of its
  "local energy" channel.
- Golden files were **re-labelled, not re-measured**: the new files were derived from the Phase 1 files
  by renaming the labels and dropping the columns each system no longer reports, and the live build then
  reproduced them bit-identically. Every dropped column (`flux_cos` everywhere, and the per-system
  always-zero `helicity_modulus` / `polyakov_loop` / `n_defects_*`) was exactly `0` in the baseline, so
  no measured number changed.
- `QXYSquare` keeps a `nPlaqs` member of its own, because its vortex loops used the base class's
  `n_plaq_variables` as their bound.

### Phase 3: RNG ownership (the one step with a re-baseline)
- Systems stop owning RNGs. `McMachine` owns one seeded `std::mt19937` and passes it to `system.randomize(rng)`
  and `system.overrelax(i, rng)`, which makes a run reproducible from a single seed.
- That changes the random stream, so verification here is **statistical**: a longer seeded run at a few T per
  system, with mean energy and observables agreeing within error bars. Then new golden files.
- Note: the same seed gives different streams on MSVC and GCC, because `std::uniform_*_distribution` is
  implementation-defined. Runs are reproducible per platform.

**As built:**
- `System` gains `randomize(std::mt19937&)`, defaulting to a no-op, and `overrelax(int, std::mt19937&)`.
  The four systems that randomized their fields now do it in `randomize`, with the same distribution as
  before (uniform(-1,1) for the two gauge systems and Spiderweb, the default [0,1) for QXYSquare).
  `XYSquare` overrides nothing: it has always started cold, and the no-op default keeps that.
- `McMachine`'s constructor seeds the generator and calls `system.randomize(rng)` straight away, so the
  initial configuration is the first thing drawn from the stream and one seed fixes the entire run.
  The systems' `seed` constructor parameters and their `rng`/`overrelax_dst` members are gone.
- `mcf_run --stats [--seeds N]` is the statistical check the phase needed, and it stays in the CLI as the
  tool for any future change to the random stream. It runs N independent chains, cooled through the same
  temperatures as the regression report, thermalizing with the step size adapting and measuring with it
  frozen. The error bar is taken **across** chains, so within-chain correlation is already accounted for.
- Result, 64 chains per system, Phase 2 build against Phase 3 build: 30 comparisons, largest deviation
  1.9 sigma, none above 3. `XYSquare` came out bit-identical, which is the expected consequence of its
  `randomize` drawing nothing. New golden files for the other three.

### Phase 4: Headless CLI for Slurm
- `SystemRegistry` plus `mcf_run`: a key=value config file with `--key=value` overrides for the system, sizes,
  couplings, `NumericalParams`, `seed` and `out`. Config parsing lives in `cli/main.cpp`, with no new dependency.
- `scripts/slurm_array.sh` as an example, and a Linux build check if you can run
  `cmake --preset linux-release && cmake --build --preset linux-release` on the cluster (I can't reach it).

**As built:**
- `SystemRegistry.h` is header-only: one `SystemEntry` per system, holding its name, its constructor
  parameters (name, default, integer-or-not) and a factory taking the values as doubles. `findSystem`
  looks one up by name. Adding a system is one entry here, and neither the CLI nor the GUI changes.
- The CLI names no system class any more: even the regression and stats cases are built through the
  registry, with their sizes spelled out so the report cannot move when a registry default does. The
  regression stayed bit-identical across that change.
- Settings are a flat key=value list: the config file first, `--key=value` overrides appended after,
  last one wins. Keys are `system`, that system's parameters, every `NumericalParams` field, and
  `seed` / `out` / `git_hash`. An unrecognised key is an **error**, not a silent no-op, because a typo
  in a config would otherwise waste a whole array job.
- `mcf_run --list` prints the registered systems and their parameters, which is also what the error
  messages point at.
- `McMachine::addProvenance(key, value)` adds `#` lines above the existing ones. The CLI uses it for
  the system name, each constructor parameter and the git hash. **The hash is passed in, not compiled
  in** (`--git_hash=$(git rev-parse --short HEAD)`, as `scripts/slurm_array.sh` does), which is what
  lets Phase 2's deferral resolve without a generated header in either build system.
- `scripts/slurm_array.sh` varies one setting per array task and writes one log per task;
  `config/spiderweb.cfg` is the example config it reads.
- Verified: both builds; regression bit-identical from both; a smoke run of **all five** systems
  producing well-formed TSVs, QXYSquare included.
- **Not verified: the Linux build.** There is no GCC on this machine and I cannot reach the cluster,
  so `cmake --preset linux-release && cmake --build --preset linux-release` is still yours to run.
  I kept to portable C++20 and added the includes MSVC supplies implicitly, but that is an argument,
  not a test.

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

**Helicity modulus — no missing `π` (checked 2026-09-22).** With `φ = πθ` the physical angle, the energy is
`-Σ cos(Δφ)` and `getSinSqrX` is `(Σ_x sin Δφ)²`, so `Γ = -E/(2N) - ⟨(Σ sin)²⟩/(T N)` is exactly the standard
`Υ = ⟨Σ_x cos Δφ⟩/N - ⟨(Σ_x sin Δφ)²⟩/(T N)` (the `/2` averages the x and y bond sums, each over N bonds).
The twist that defines Υ is a phase on the spin, `φ → φ + δ`, not a shift of `θ`, so the derivative brings down no
`π` and both terms share one normalisation. The earlier entry in this appendix was wrong and has been removed.

### Still open, for you (I will not touch these)
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
