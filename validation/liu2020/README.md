# Liu et al. (2020): carbamazepine cooling crystallization

> Y. C. Liu, D. Acevedo, X. Yang, S. Naimi, W.-L. Wu, N. Pavurala, Z. K. Nagy,
> T. F. O'Connor, *Population Balance Model Development Verification and
> Validation of Cooling Crystallization of Carbamazepine*, Cryst. Growth Des.
> 2020, **20**, 5235-5250.

Seeded batch cooling of carbamazepine (form III) in ethanol, 350 mL, 350 rpm.
The four validation experiments Val1-Val4 of the paper, run with the 0-D
solver `crystalBatchPbeFoam` and compared with the measured concentration.

## Cases

Every case is a copy of `examples/liu2020_batch` (same kinetics: growth GM5,
secondary nucleation BM5, dissolution DM1, Table 10 of the paper). Only the
seed, the temperature programme and the end time differ:

| case | seed | sieve cut | temperature programme | end |
|---|---|---|---|---|
| `val1` | 1.28 g | 75-125 um | 35 -> 10 degC at 0.2 degC/min, hold | 150 min |
| `val2` | 0.96 g | 75-125 um | 35 -> 10 degC at 0.4 degC/min, hold | 90 min |
| `val3` | 1.60 g | 75-212 um | 35 -> 10 degC at 0.2 degC/min, hold | 150 min |
| `val4` | 1.60 g | 75-212 um | 35 -> 10 degC at 0.2 degC/min, then 10 -> 45 degC at 1 degC/min | 170 min |

All start from C = 40 mg/mL. The seed has a uniform number density over the
sieve cut, scaled so that `kv rhoCrystal m3` equals the seed mass per volume
(seed mass / 350 mL); this gives the moments in `0/moment.*`.

Units as in the example: L in um, C in mg/mL, time in s.

`refdata/<case>.csv` holds the curves digitized from the paper: columns 1-2
are the experiment, columns 3-4 the paper's model (time in min, C in mg/mL).

## Run

With the OpenFOAM environment sourced and the solver built:

```bash
./Allrun              # all four cases (a few seconds each)
./Allrun val4         # one case
./Allclean            # remove all run output
```

## Analyze

```bash
python3 analyze.py            # every case that has been run
python3 analyze.py val1 val4  # only these
```

It prints, per case, the RMSE of the concentration against the experiment for
OpenFOAM and for the paper's model, plus end values (C, m0, d10, d43), and
writes:

* `results/liu2020.png`: concentration (experiment, paper model, OpenFOAM,
  solubility) and the mean sizes d10, d43 over time, one row per case;
* `results/summary.csv`: the printed table.

Run output (`postProcessing/`, time directories, logs, mesh) is ignored by git.
