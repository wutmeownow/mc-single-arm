# 3He SF fit tables

Place the CSV tables required by `GETSF_F1F2fit` in this directory:

- `Table_3He_F1F2_SF1.csv`
- `Table_3He_F1F2_SF2.csv`
- `Table_3He_F1F2_SF3.csv`
- `Table_3He_F1F2_SF4.csv`
- `Table_3He_F1F2_SF5.csv`

The code now probes both:

1. `sf_tables/` (runtime working directory)
2. `src/interp/sf_tables/` (repo-local fallback)

If none are found, `GETSF_F1F2fit` returns `STAT=.false.` and caller fallback logic is used.


## Automatic extraction during build

When you run `make` in `src/`, the `prepare_sf_tables` step now checks whether
the CSVs already exist. If not, it attempts to extract one of:

- `interp/F1F221_3He_XZ_20250828_tables.tar.gz`
- `interp/F1F221_3He_XZ_20250828_tables.tar`
- `interp/F1F221_3He_XZ_20250828_tables.tar.*` (split tar parts)

into `interp/`, so `interp/sf_tables/` is populated automatically when the
archives are present.