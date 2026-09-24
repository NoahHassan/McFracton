# Spiderweb specific heat: derivative vs. fluctuations

**Question.** Do the two routes to the specific heat agree over a range of temperatures?

- **derivative route** `C = d⟨H⟩/dT`, from the logged `Energy` column
- **fluctuation route** `C = (⟨H²⟩ − ⟨H⟩²)/T²`, from the logged `DEnergy` column

**Answer.** They agree, at every size and coupling tested, over three decades of temperature.
There is no bug in the Spiderweb physics or in the Monte Carlo.

**But there is a real trap, and it produces exactly the symptom you described.** The `DEnergy`
column no longer means what it meant before the refactor. Using the formula that was correct for
your older log files gives a specific heat that **rises as `1/T` all the way down instead of
saturating**. Details and the fix are in section 4.

Data, scripts and the figure: `results/specific_heat_check/` (gitignored).
Figure: `results/specific_heat_check/specific_heat_check.png`.

---

## 1. What was measured

Annealing runs through the usual geometric schedule (`t_max = 100`, `t_fac = 0.9`, down to
`T = 0.05`), with the step size adapting during thermalisation and frozen during measurement, as
the code already does.

| Run | L, N_t | KU | seeds | `n_measurements` |
|---|---|---|---|---|
| `L4_s1..s4` | 4, 4 | 0.5 | 4 | 100 |
| `L4_long` | 4, 4 | 0.5 | 1 | 300 |
| `L8_n100_s1..s2` | 8, 8 | 0.5 | 2 | 100 |
| `L4_KU0.1` | 4, 4 | 0.1 | 1 | 100 |
| `L4_KU2.0` | 4, 4 | 2.0 | 1 | 100 |
| `L16_default` | 16, 16 | 0.5 | 1 | 10 (config default) |

`d⟨H⟩/dT` is a central difference on the (non-uniform) temperature grid. The variance of a single
measurement is recovered from the error column as `Var = n_measurements · DEnergy²` — see section 4.

## 2. Result: the two routes agree

Ratio `(d⟨H⟩/dT) ÷ (Var(H)/T²)` over `T ∈ [0.05, 3]`:

| Run | mean | median | min | max |
|---|---|---|---|---|
| L=4, KU=0.5, 4 seeds | 1.023 | 1.010 | 0.796 | 1.249 |
| L=4, KU=0.5, 3× statistics | 1.056 | 1.039 | 0.886 | 1.425 |
| L=8, KU=0.5, 2 seeds | 1.038 | 1.026 | 0.774 | 1.434 |
| L=4, KU=0.1 | 1.031 | 0.985 | 0.758 | 1.584 |
| L=4, KU=2.0 | 1.081 | 1.076 | 0.680 | 1.696 |
| L=16, KU=0.5, default config | 1.476 | 1.205 | 0.444 | 5.577 |

The scatter is statistical: it shrinks when the statistics grow, and it is centred on 1. The last
row is the exception and it is not a physics problem — see section 5.1.

Both curves have the **same shape**: a peak near `T ≈ 0.3` and a flat plateau below it. Neither
diverges, neither saturates while the other diverges.

## 3. Three independent checks that there is no bug

### 3.1 The proposal is consistent with the energy

If `proposeUpdate` did not return the exact change in `getEnergy()`, the chain would sample
`exp(−H′/T)` for some `H′ ≠ H`, and `d⟨H⟩/dT = Var(H)/T²` would fail *by construction*. This is
precisely the bug you fixed in `AbelianGaugeCube::getLocalEnergy_x/y` in commit `94bdb49`.

`results/specific_heat_check/probe.cpp` compares `proposeUpdate(i, δ)` against
`getEnergy()` before and after `applyUpdate(i, δ)`, for 4000 random full-size proposals per size:

```
L=4   Nt=4   KU=0.5    worst abs err 1.421e-14   worst rel err 9.024e-12
L=6   Nt=6   KU=0.5    worst abs err 1.599e-14   worst rel err 4.653e-12
L=6   Nt=4   KU=0.5    worst abs err 1.354e-14   worst rel err 2.574e-11
L=7   Nt=5   KU=0.5    worst abs err 3.553e-14   worst rel err 7.279e-11
L=8   Nt=6   KU=0.5    worst abs err 3.730e-14   worst rel err 2.223e-11
L=8   Nt=6   KU=0.1    worst abs err 1.501e-13   worst rel err 3.097e-11
L=8   Nt=6   KU=2      worst abs err 9.770e-15   worst rel err 8.212e-10
L=5   Nt=3   KU=0.5    worst abs err 9.603e-15   worst rel err 2.122e-11
```

Machine precision throughout. The odd sizes matter: the Spiderweb stencil reaches ±2 sites in x and
y, so on `L = 4` the site `r+2x` *is* the site `r−2x` and a swapped-sign stencil error would cancel
silently. `L = 5, 6, 7, 8` with `L ≠ N_t` close that loophole.

### 3.2 The historical bug was real, and it is gone

At `94bdb49` — the commit whose message says spiderweb "still has a specific heat bug" —
`Spiderweb::proposeSiteFlip` computed

```cpp
auto e_terms_xx = getElectricTerms_xx(index);
auto e_terms_xy = getElectricTerms_xy(index);
```

i.e. only the terms **anchored at the changed site**. A variable also enters the terms anchored at
`r−x`, `r−y`, `r−2x`, `r−2y`, `r−t` … , which were all missing, so `ΔE ≠ ΔH`. That is the same
class of bug as the Cube one, in the same commit. You rewrote it into the present multi-anchor
form (the `accumulate_xx` / `accumulate_xy` / `accumulate_b` lambdas) before `fa44d03`, and the
refactor carried that through bit-identically. **The bug your commit message refers to has already
been fixed.**

### 3.3 The low-temperature plateau matches equipartition, with no Monte Carlo involved

`results/specific_heat_check/hessian.cpp` quenches the system to a local minimum and builds the
Hessian of `H` there (from `proposeUpdate`, which section 3.1 showed is exact). The spectrum has a
clean gap — the mode count is identical whether the threshold is `|λ| > 0.1` or `|λ| > 1`, and
identical for independent quenches:

| | variables | massive modes | flat directions |
|---|---|---|---|
| L=4, N_t=4 (seed 1) | 192 | 126 | 66 |
| L=4, N_t=4 (seed 2) | 192 | 126 | 66 |
| L=8, N_t=8 (seed 3) | 1536 | 1022 | 514 |

The flat directions count exactly `N_sites + 2` (66 = 64+2, 514 = 512+2) — the rank-2 gauge
redundancy, one gauge function per site, plus two global modes. Classical equipartition then gives

> **C(T → 0) = (3N_sites − N_sites − 2)/2 = N_sites − 1**

with no free parameter. Measured plateau (mean of the five coldest points):

| Run | `d⟨H⟩/dT` | `Var(H)/T²` | `N_sites − 1` |
|---|---|---|---|
| L=4, KU=0.5, 4 seeds | 64.9 | 65.7 | 63 |
| L=4, KU=0.5, 3× stats | 69.0 | 68.3 | 63 |
| L=8, KU=0.5, 2 seeds | 538.0 | 514.9 | 511 |
| L=4, KU=0.1 | 67.6 | 58.6 | 63 |
| L=4, KU=2.0 | 66.9 | 71.0 | 63 |

Both routes land on the prediction, and the plateau is independent of `KU` — as it must be, since
`KU` changes the curvatures but not the number of massive modes.

## 4. The real trap: `DEnergy` is no longer a variance

**This is the one thing that reproduces your symptom.**

Before the refactor (`fa44d03` and earlier), `McMachine::Measure` logged

```cpp
s_sqr.energy += std::pow((observables_T[n].energy - means.energy), 2.0f) / N;   // = s² = (1/N)Σ(E−Ē)²
```

so the `DE` column **was** the variance of a single measurement, and

```
C = DE / T²
```

was correct for those files.

Phase 0 changed the error columns to the standard error of the mean (`REFACTOR_PLAN.md`,
Appendix A.4 — an agreed change). `McMachine::Measure` now logs

```cpp
return std::sqrt(sum_sqr / (double(N) * double(N - 1)));   // = s/√N
```

The same formula applied to a current log therefore evaluates to

```
DEnergy / T²  =  s / (√N T²)  =  (T√C) / (√N T²)  =  √C / (√N · T)
```

Since `C → const` as `T → 0`, this **grows like `1/T` forever** while the true `C` flattens. On the
L=4 data, `DEnergy/T²` climbs monotonically from 0.76 at `T = 0.97` to 17.0 at `T = 0.05`, with no
peak and no plateau, and it follows the `1/T` guide line in the bottom row of the figure. That is
"one diverges at low T while the other saturates", produced entirely by the column's change of
meaning.

### The correct reading today

```
Var(H)   = n_measurements · DEnergy²
C        = n_measurements · DEnergy² / T²
```

`n_measurements` is in the `#` header of every log file, so it can be parsed rather than hard-coded.

### Two related snags

- **Sign.** Your message writes the fluctuation formula as `1/T²(⟨H⟩² − ⟨H²⟩)`. That is *minus*
  the variance. If a script implements it literally the curve comes out negative everywhere and
  vanishes from a log plot.
- **Column rename.** The value columns kept their old spellings, but the error columns went from
  `DE`/`DHM`/`Dna` to `DEnergy`/`DHelicity Modulus`/… A script that looks the column up **by name**
  now raises a `KeyError`, which is safe. A script that reads it **by position** silently gets the
  new quantity — and the position moved too, because each system now logs only the columns it
  actually computes (Spiderweb logs `T, Energy, DEnergy, autocorrelation, n_sweeps, acceptance`,
  where the old format always had all five observables).

### Suggested fix

My preference, in order:

1. **Say what the column is, in the file.** Add one provenance line in `McMachine::StartSimulation`,
   next to the existing `# variables` line:

   ```cpp
   logfile << "# error_columns\tstandard_error_of_mean\tn_measurements\t" << params.n_measurements << '\n';
   ```

   This costs nothing, breaks no parser that reads with `comment='#'`, and does not touch the
   regression goldens (`--regression` drives `Sweep`/`measure` directly and never opens an annealing
   log). It makes every future file self-describing, including the ones already on the cluster.

2. **Log the variance as well**, as a second derived column, so the specific heat needs no
   reconstruction at all. This *does* change the log format and the column count, so under rule 4 it
   needs your agreement first, and it would want new goldens if it ever touched the regression path.

3. **Leave the code alone** and fix it in the analysis script. Fine if your analysis lives in one
   place, but the old `.txt` files and the new ones then need different formulas with nothing in
   either file to say which is which — which is how this happened.

I have not implemented any of these. Option 1 is a two-line change I can make whenever you say so.

## 5. Other things worth knowing (not bugs, but they cost me time)

### 5.1 `n_measurements = 10` is too few for a fluctuation-based specific heat

The variance of `N` samples has a relative uncertainty of about `√(2/(N−1))` — 47 % at the config
default of `N = 10`. That is the whole story of the L=16 row in section 2 (ratio range 0.44–5.58):
the derivative route is smooth because the annealing curve is smooth, while the fluctuation route
scatters wildly. Nothing is wrong with it, there is just no information in 10 samples.

For a usable `C` from fluctuations, `n_measurements` wants to be ~100 or more. At `N = 100` the
scatter dropped to 0.8–1.25, and at `N = 300` to 0.89–1.43 on a single seed.

### 5.2 The first temperature can cost more than the whole rest of the run

`McMachine`'s constructor initialises

```cpp
current_measurement_sweeps(params.initial_therm_sweeps)
```

ignoring `max_measure_sweeps`, which only takes effect after the *first* `Measure` has run and the
autocorrelation has been evaluated. So the first temperature performs
`n_measurements × initial_therm_sweeps` sweeps regardless of the `max_measure_sweeps` in the config.
With `n_measurements = 300` and `initial_therm_sweeps = 500` that is 150 000 sweeps for `T = 100`,
where the system is trivially disordered — it dominated the wall time of two of my runs before I
noticed.

Suggestion (pure performance, no numbers change unless a run actually hits the clamp):

```cpp
current_measurement_sweeps(std::min(params.initial_therm_sweeps, params.max_measure_sweeps))
```

This *would* change logged numbers for any run where `initial_therm_sweeps > max_measure_sweeps`,
so it is your call, not mine.

### 5.3 `updates_per_sweep` is an absolute count, not per-variable

The default 5000 is smaller than the variable count as soon as the lattice is reasonably large:
`L = 16, N_t = 16` has 12 288 variables, so one "sweep" touches about 40 % of them. "20 sweeps
between measurements" then means 8 lattice sweeps, and the autocorrelation time in the log is
measured in units that change meaning with system size. Not wrong, but worth remembering when
comparing `n_sweeps` across sizes, and worth a line in `config/spiderweb.cfg`.

## 6. Reproducing this

```sh
# the scans (results/ is gitignored)
build/windows-release/Release/mcf_run.exe --config config/spiderweb.cfg \
    --linear_size=4 --temporal_size=4 --t_min=0.05 --updates_per_sweep=1000 \
    --n_measurements=100 --max_measure_sweeps=20 --max_therm_sweeps=200 \
    --initial_therm_sweeps=50 --seed=1 --out=results/specific_heat_check/L4_s1.txt

# the comparison table and the figure
python results/specific_heat_check/analyse.py "results/specific_heat_check/L4_s*.txt"
python results/specific_heat_check/plot.py results/specific_heat_check
```

`probe.cpp` and `hessian.cpp` are standalone; they link against `mcf_core.lib` and are deliberately
**not** part of the solution or of `CMakeLists.txt`:

```sh
cl /nologo /std:c++20 /O2 /MD /EHsc /IMcFracton\core\include \
   results\specific_heat_check\probe.cpp build\windows-release\Release\mcf_core.lib
```
