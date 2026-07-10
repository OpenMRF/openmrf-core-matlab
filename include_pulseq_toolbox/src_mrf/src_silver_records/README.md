# Silver-dictionary records

A small toolkit for producing **silver records**: self-describing, interoperable
bundles of an ideal IR-FISP (+ optional T2-prep) EPG signal dictionary together
with everything needed to reproduce and cross-validate it. A "silver" record is
a reproducible reference standard used to gauge the *commensurability* of a
sibling simulator (e.g. a Python EPG port) against the MATLAB reference.

Signal generation is delegated entirely to `MRF_sim_irfisp_epg`; this toolkit
only packages, serializes and compares.

## What a record contains

One self-contained JSON file (schema `openmrf.silver/1`):

| Section       | Contents |
|---------------|----------|
| `manifest`    | array shapes, **dim order**, units, **complex & RF-handedness conventions**, `M0` normalisation, demodulation note |
| `inputs`      | `FAs` (real+imag), `TRs`, `TE`, prep spec, spoiler/phase options, `P.{T1,T2,...}`, tissue mode + labels |
| `outputs`     | `Mxy_real`, `Mxy_imag`  `(NR x N_dict)` |
| `provenance`  | timestamp, git hash + dirty flag, machine id, host/user, OS, MATLAB version, cores |
| `integrity`   | SHA-256 of the numeric payload (metadata excluded, so identical numbers -> identical hash) |

## Functions

| File | Purpose |
|------|---------|
| `make_silver_record(P, FAs, TRs, opt, meta)` | run the sim, assemble the record struct (no I/O) |
| `write_silver_record(rec, path)`             | serialize to `.json` |
| `read_silver_record(path[, verify])`         | load, reshape, reconstruct complex, verify checksum |
| `compare_silver_records(A, B[, tol])`        | commensurability gauge: structure check + per-atom RMSE / max-abs / correlation |
| `default_silver_tissues([names])`            | curated nominal (T1,T2) tissue table |
| `silver_checksum`, `silver_sha256`           | payload hashing helpers |
| `silver_preset(name)`                        | resolve a named record-flavour preset -> build spec |
| `list_silver_presets()`                      | enumerate available presets (metadata only) |
| `make_silver_record_from_preset(name[, dir])`| resolve preset -> build -> write, in one call |
| `demo_silver_record.m`                       | driver: build + write + round-trip a record |
| `demo_silver_presets.m`                      | driver: list presets + build the default by name |

## Usage

```matlab
install_OpenMRF          % once, to put everything on the path
demo_silver_record       % builds a record into <repo>/silver_records/
```

Generated records land in `<repo>/silver_records/` (git-ignored); regenerate
them from the committed driver rather than committing the artifacts.

## Presets (record flavours)

A **preset** is a named, deterministic recipe for a record's *inputs* (pattern,
tissues, sequence options) recalled by a memorable `adj-adj-noun` codename or a
short alias. Definitions live as language-neutral JSON in `presets/`; the
pattern is referenced by name (`"yun"`) and resolved via `MRF_get_FAs_TRs`. Each
produced record is stamped with `rec.preset = {name, title, spec_sha256}`, where
`spec_sha256` fingerprints the numeric excitation recipe (FAs, TRs, TE, T1, T2)
so a record self-identifies which flavour built it.

| Codename | Alias | Flavour |
|----------|-------|---------|
| `calm-silver-heron`   | `default` | IR-FISP, canonical Yun 1000-TR, 10 named tissues (baseline) |
| `brisk-teal-marlin`   |           | IR + mid-train T2-prep (tau=60 ms @ readout 500), Yun 1000-TR, named tissues |
| `dense-amber-lattice` |           | IR-FISP, Yun 1000-TR, dense log-grid (T2<T1) tissues |

```matlab
list_silver_presets();                              % see what's available
make_silver_record_from_preset('default');          % build + write by name
spec = silver_preset('brisk-teal-marlin');          % or resolve, then build manually
```

To add a flavour, drop a new `presets/<codename>.json` (schema
`openmrf.silver.preset/1`); no code change needed. A sibling tool can read the
same file to resolve the same codename to the same configuration.

## Cross-tool validation

Have the sibling tool emit the identical JSON structure for the same inputs, then:

```matlab
rep = compare_silver_records('silver_records/reference.json', 'other_tool.json', 1e-4);
```

The report gates on schema / layout / units / conventions first, then aligns
tissue atoms by their `(T1,T2)` key and quantifies signal agreement.
