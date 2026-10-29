# sample-input

This app reads the output of `processor-fcs-prep` (and, in the workflow, of `processor-peacoqc`,
`processor-cytonorm` or `processor-cycombine`): `fcs/<sample_id>.fcs` with arcsinh-transformed
marker channels and `metadata/{samples.csv, panel.csv, prep.json}`.

## Synthetic cohort (no download)

```sh
make fixture PYTHON=.venv/bin/python   # data/fixture/{linked,flat}/
make run PYTHON=.venv/bin/python       # Docker over data/fixture/linked/ -> data/output/raw/
```

4 batches × {unstim, stim} × 2 samples + one reference sample (`control=ref`) per batch,
3000 events each, 5 markers. Batch `b3` is shifted by +0.6 on M1 and M3 (flagged), stimulation
raises M2 in every batch (biology, not flagged), and prep.json says batch `b4` lacked M6.

## Real data: the CytoBatchFlagR spectral flow dataset

Koenig Z, Andorf S, Tilburgs T, Eswar S. CytoBatchFlagR: Spectral flow cytometry dataset.
Zenodo 2025, https://doi.org/10.5281/zenodo.15388817 (CC BY 4.0): 24 FCS files, 6 batches ×
{tcu1, tcu2 unstimulated; tcs3, tcs4 stimulated}, 1.46 GB.

1. Download the FCS files from Zenodo (each file's link is
   `https://zenodo.org/api/records/15388817/files/<name>/content`).
2. Run `processor-fcs-prep` on them with the sample sheet (`file, sample_id, batch, condition,
   control`) to get `fcs/` and `metadata/`.
3. Run this app with `REFERENCE_COLUMN=control REFERENCE_VALUES=tcu1`, or make a smaller copy
   with a planted batch effect first:

```sh
make fixture-demo PYTHON=.venv/bin/python PREP=/path/to/fcs-prep/output
# -> data/fixture-demo/: 20k events per file, batch03 shifted +0.4 on CD3, CD4, CD8
```
