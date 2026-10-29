# processor-cytometry-batch-assess

How much do the batches of a cytometry cohort differ, and on which markers? This app measures
batch effects in a flow or mass cytometry cohort so a researcher can see **which batches and
markers drift** and, run again after a correction, **whether the correction helped**.

For every batch and marker it computes the **Earth Mover's Distance (EMD)** between expression
distributions: the amount of "intensity" each cell would have to move, on average, to turn one
distribution into the other (in the same units as the transformed data). EMD was introduced for
comparing marker expression in cell populations by Orlova et al. (2016); cyCombine (Pedersen et
al. 2022) uses the reduction in EMD between batches to evaluate batch correction, and
CytoBatchFlagR (Eswar et al. 2026), made by the authors of our demo dataset, uses EMD among its
batch flags. This app is an independent Python implementation of the same idea (no code is
taken from those tools; CytoBatchFlagR is GPL-3.0).

It is a step of the **Multi-batch Harmonization & QC** workflow:
`fcs-prep → peacoqc → {batch-assess (LABEL=raw), cytonorm → batch-assess (LABEL=cytonorm),
cycombine → batch-assess (LABEL=cycombine)} → cytometry-report`. All three runs link into one
input folder of the report, so **everything is written under `OUTPUT_DIR/<LABEL>/`**.

## Method

1. **Per sample, per marker: a distribution.** Each FCS file is read once (one at a time),
   up to `CELLS_PER_SAMPLE` cells are drawn without replacement (seeded per sample from `SEED`
   and the sample id, so a sample's draw doesn't depend on the rest of the run), and each
   marker's transformed values are binned on a fixed grid (−10 to 20, bins of 0.005). Values
   are used as they are: fcs-prep has arcsinh-transformed the markers, and downstream steps
   never re-transform.
2. **Pooling.** A batch's distribution is the equal-weight mixture of its samples; "all batches"
   is the equal-weight mixture of the batches, so a batch with more or bigger files doesn't
   dominate.
3. **EMD.** `scipy.stats.wasserstein_distance` (1-D Wasserstein-1 = EMD) between the binned
   distributions, on the bin centres weighted by bin mass. Binning moves each value by at most
   half a bin, so the result is within 0.005 of the EMD on the raw values (tested).
4. **Two comparisons.**
   - **Reference design** (`REFERENCE_COLUMN` / `REFERENCE_VALUES` set): reference samples are
     aliquots of the same sample run in every batch (in the demo, `control=tcu1`). Each batch's
     reference sample is compared with the reference samples of all batches pooled. The cells
     are the same everywhere, so any distance is **technical variation**. With several
     reference values, each is compared with its own pool.
   - **Batch vs batch** (always): within each group of `CONDITION_COLUMN`, a batch's cells are
     compared with the cells of all *other* batches pooled. Stratifying by condition keeps a
     biological difference (stimulated vs unstimulated) from being counted as batch effect, as
     long as the conditions are spread over the batches.
5. **Flags.**
   - A batch × marker distance is flagged when it is above **median + `MAD_K` × MAD** of the
     same marker's distances across batches (per reference value, or per condition; MAD scaled by
     1.4826 so it estimates the SD for normal data; at least 3 batches needed) **and** above
     `MIN_EMD`.
   - **Missing markers**: fcs-prep keeps only channels present, with the same marker, in every
     sample and lists in `prep.json` what each sample had that was dropped. A dropped marker is
     flagged for the batches that lack it (when at most half the batches lack it), or for the
     batches that have it (◆, "panel differs", when most batches lack it). A panel marker absent
     from an FCS file is flagged the same way. In the demo, batches 05–06 have CCR4 on APC-A
     instead of GITR and CD28 on PerCP-A instead of CD163.
   - A **batch** is flagged when at least `BATCH_FLAG_FRACTION` of the markers are flagged for it.
   - **Per-batch score**: the mean EMD over markers, from the reference comparison when there is
     one, else from the batch-vs-batch comparison. Lower is better; compare it across LABELs.

## Ports

| Port | Kind | Files | Notes |
|---|---|---|---|
| `fcs` | input (subfolder) | `<sample_id>.fcs` | FCS 3.1 float32, marker channels arcsinh-transformed (fcs-prep, peacoqc, cytonorm or cycombine output). Other files in the folder are ignored. |
| `metadata` | input (subfolder) | `samples.csv`, `panel.csv`, `prep.json` | fcs-prep's metadata. `samples.csv` needs `sample_id`, `batch`, and the `CONDITION_COLUMN` / `REFERENCE_COLUMN` columns used; `panel.csv` rows with `type=marker` are assessed; `prep.json` is optional. |
| `assessment` | output | `<LABEL>/…` (below) | The only output, so it is the whole of `OUTPUT_DIR`; the app writes into `OUTPUT_DIR/<LABEL>/` so several runs can share a folder downstream. |

Without links (an older compute node, or files dropped together) the same files are found by
name anywhere under `INPUT_DIR` (symlinks followed). Without `samples.csv` the batch is each file's
`$DATE` keyword and there is one condition; without `panel.csv`, channels with a `$PnS` label
(other than scatter/time) are assessed.

## Parameters

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `LABEL` | string | `raw` | Name of this assessment, `[a-z0-9-]+` (e.g. `raw`, `cytonorm`, `cycombine`); the output folder. Anything else fails the run. |
| `REFERENCE_COLUMN` | string | (empty) | `samples.csv` column marking reference samples (e.g. `control`). Empty: no reference comparison. |
| `REFERENCE_VALUES` | string | (empty) | Comma-separated values of that column that are references (e.g. `tcu1`). Set both or neither. |
| `CONDITION_COLUMN` | string | `condition` | Column whose groups are compared separately in batch vs batch. |
| `CELLS_PER_SAMPLE` | integer | `50000` | Cells drawn per sample (seeded); `0` = all. |
| `MAD_K` | number | `3` | K in median + K × MAD. |
| `MIN_EMD` | number | `0.1` | Distances at or below this are never flagged. On an arcsinh scale 0.1 is roughly a 10 % intensity shift of the positive cells. |
| `BATCH_FLAG_FRACTION` | number | `0.2` | A batch is flagged when at least this fraction of the markers is flagged for it. |
| `MARKERS` | string | (empty) | Comma-separated markers or channels to assess; empty = every panel marker. |
| `MAX_PLOTS` | integer | `6` | Density plots for at most this many markers. |
| `SEED` | integer | `42` | Subsampling seed. |

Demo settings: `REFERENCE_COLUMN=control`, `REFERENCE_VALUES=tcu1`.

## Outputs (`OUTPUT_DIR/<LABEL>/`)

| File | Contents |
|---|---|
| `emd_reference.csv` | `batch, marker, emd, reference, threshold, flagged, n_samples, n_cells`: each batch's reference sample vs the pooled references. Header only without a reference design. |
| `emd_batches.csv` | `batch, condition, marker, emd, threshold, flagged, n_samples, n_cells`: each batch vs the other batches, within a condition. |
| `flags.csv` | `level` (`marker` / `batch`), `comparison` (`reference`, `batches`, `missing`, `summary`), `batch, condition, reference, marker, emd, threshold, reason`. |
| `summary.json` | `label`, `n_batches`, `n_markers`, `n_samples`, `median_emd` (+ `median_emd_source`), per-comparison stats (`n`, `median_emd`, `mean_emd`, `max_emd`, `n_flagged`), `per_batch` (`score`, `mean_emd_reference`, `mean_emd_batches`, `n_flagged_markers`, `flagged_markers`, `flagged`), `flagged_batches`, `flags`, `markers`, `parameters`, `figures`. |
| `heatmap_reference.png` | Batch × marker EMD of the reference samples; flagged cells outlined in red with their value, ✕ missing marker, ◆ marker only in these batches. |
| `heatmap_batches.png` | The same for batch vs batch, one panel per condition. Both heatmaps of a run share one colour scale. |
| `density_<marker>.png` | Ridge plots of the most-flagged markers (then the largest distances): one row per batch (the reference sample, or each condition), the pooled distribution dashed, flagged batches marked ▲ in red. |
| `provenance.json` | App and tool versions, parameters, method constants, inputs seen, per-sample status (events, cells used, markers absent, values outside the grid, seconds), timings, peak RSS. |

A sample that cannot be read is logged and skipped; the run fails (exit 1) only when no sample
could be read or none was found, and exits 2 on bad parameters.

## Results on the demo data

fcs-prep output of the CytoBatchFlagR dataset (24 files, 6 batches, 10.4 M events, 25 markers
kept), `REFERENCE_COLUMN=control REFERENCE_VALUES=tcu1`, defaults otherwise (`LABEL=raw`):
median reference EMD 0.135 (batch vs batch 0.155); per-batch scores 0.29 (batch01), 0.10, 0.14,
0.13, 0.14, 0.20 (batch06). Flagged batches: **batch01** (10 markers, e.g. CD45RA 0.99, CD3 0.58,
CD27 0.52 in the reference sample), **batch06** (16, incl. GZMB, CD8, CD14) and **batch05**
(5, mostly the panel change: GITR, CD163 missing, CCR4 only there). With a planted +0.4 shift
on CD3, CD4 and CD8 in batch03 (`make fixture-demo`, 20 k events per file), batch03's
reference EMDs rose from 0.22 / 0.20 / 0.13 (same subsample, no shift) to 0.51 / 0.52 / 0.46,
and CD4 and CD8 were flagged; CD3 was not (threshold 0.89), because batch01 is itself far off on
CD3 and two outliers in six batches inflate the MAD.

## Runtime

`computeTypes: [standard, lambda]`, 1 vCPU / 2 GiB, 15 min timeout. Memory doesn't grow with
cohort size beyond the histograms (~1 MB per sample at 25 markers): files are read one at a time
and each is reduced to histograms before the next is read.

Measured in the linux/amd64 image under emulation on an M-series Mac (Docker, `--memory 1g`),
on the full demo cohort (fcs-prep output, 24 files, 1.4 GB, 10.4 M events), reference design:
**peak RSS 346–353 MB, wall time 5.2–5.6 s** over runs (read 2.0 s, compare 0.1 s, figures 1.5–2 s). The
largest file (603 k events × 33 channels) dominates memory: FlowIO holds a whole file. The
cgroup's `memory.peak` reaches the limit only through page cache of the files read
(reclaimable; `anon` stayed at ~5 MB at exit). Image: 168 MB.

Lambda: 2 GiB leaves 5× headroom over the measured peak, and the time scales at roughly
0.1–0.2 s per file, so cohorts of a few thousand files fit in the 15-minute limit.

## Running it locally

```sh
python3 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
make test PYTHON=.venv/bin/python            # unit tests on the synthetic cohort, no Docker
make run PYTHON=.venv/bin/python             # Docker over data/fixture/linked -> data/output/raw/
make fixture-demo PYTHON=.venv/bin/python    # subsampled demo data with a planted shift
```

By hand:

```sh
docker build --platform linux/amd64 -t processor-cytometry-batch-assess .
docker run --rm --platform linux/amd64 --memory 1g \
  -e INPUT_DIR=/data/input -e OUTPUT_DIR=/data/output \
  -e LABEL=raw -e REFERENCE_COLUMN=control -e REFERENCE_VALUES=tcu1 \
  -v "$PWD/data/fixture-demo:/data/input:ro" -v "$PWD/data/output:/data/output" \
  processor-cytometry-batch-assess
```

The same image runs as a Lambda function: with `AWS_LAMBDA_RUNTIME_API` set, `processor/entry.py`
hands over to the runtime interface client with `processor.handler.handler`, which takes
`inputDir`, `outputDir` and the parameters from the event.

## Limitations

- EMD is computed on each marker's whole distribution, one marker at a time. A shift confined
  to one cell population, or a change in the joint distribution, can go unnoticed when the
  marginal stays similar; cyCombine computes EMD per marker per SOM cluster, and CytoBatchFlagR
  adds cluster-level and positive/negative-population measures. Cluster-level assessment is not
  done here.
- The batch-vs-batch comparison assumes the conditions (and the donors) are similar across
  batches. If batch 3 holds the sick patients, biology is counted as batch effect. The reference
  comparison doesn't have this problem, which is why reference samples are recommended.
- The flag threshold is relative to the cohort's own batches: with few batches (three to six)
  and more than one bad batch, the MAD grows and real shifts can stay below the threshold (CD3
  above). With two batches nothing is flagged; read the distances.
- Distances are in transformed units, so they compare across markers and runs only when the
  transform (cofactor) is the same — as it is within one workflow.
- Values outside −10…20 are counted in the edge bins and a warning is logged per sample (it
  would mean the markers weren't transformed).

## Licences

Apache-2.0 for this app (see `LICENSE`). Bundled components and their licences: `THIRD_PARTY.md`.

## Citations

Earth Mover's Distance for comparing marker expression in cell populations (the method):

> Orlova DY, Zimmerman N, Meehan S, Meehan C, Waters J, Ghosn EE, Filatenkov A, Kolyagin GA,
> Gernez Y, Tsuda S, Moore W, Moss RB, Herzenberg LA, Walther G.
> **Earth Mover's Distance (EMD): A True Metric for Comparing Biomarker Expression Levels in Cell
> Populations.** PLoS One 11(3):e0151859 (2016). https://doi.org/10.1371/journal.pone.0151859,
> PMID 27008164.

EMD reduction as the measure of batch correction:

> Pedersen CB, Dam SH, Barnkob MB, Leipold MD, Purroy N, Rassenti LZ, Kipps TJ, Nguyen J,
> Lederer JA, Gohil SH, Wu CJ, Olsen LR. **cyCombine allows for robust integration of single-cell
> cytometry datasets within and across technologies.** Nat Commun 13(1):1698 (2022).
> https://doi.org/10.1038/s41467-022-29383-5, PMID 35361793.

The same idea, as an R tool, by the authors of the demo dataset (not used in code; GPL-3.0,
https://github.com/AndorfLab/CytoBatchFlagR):

> Eswar S, Koenig ZT, Tursi AR, Cobeña-Reyes J, Tilburgs T, Andorf S. **CytoBatchFlagR: A
> Comprehensive Framework to Objectively Assess High-Parameter Cytometry Data for Batch
> Effects.** Cytometry A 109(4):241-256 (2026). https://doi.org/10.1002/cyto.a.70024,
> PMID 41968341.

The EMD itself:

> Rubner Y, Tomasi C, Guibas LJ. **The Earth Mover's Distance as a Metric for Image Retrieval.**
> International Journal of Computer Vision 40:99-121 (2000).
> https://doi.org/10.1023/A:1026543900054 (not in PubMed; checked on Crossref).

SciPy (`wasserstein_distance`):

> Virtanen P, et al. **SciPy 1.0: fundamental algorithms for scientific computing in Python.**
> Nat Methods 17:261-272 (2020). https://doi.org/10.1038/s41592-019-0686-2, PMID 32015543.

FCS reading: FlowIO, https://github.com/whitews/FlowIO (no paper of its own).

Demo data:

> Koenig Z, Andorf S, Tilburgs T, Eswar S. **CytoBatchFlagR: Spectral flow cytometry dataset.**
> Zenodo (2025). https://doi.org/10.5281/zenodo.15388817 (CC BY 4.0; checked on DataCite).

Every PMID above was checked on PubMed (E-utilities `esummary`: title, first author, year,
journal, volume, pages and DOI match).
