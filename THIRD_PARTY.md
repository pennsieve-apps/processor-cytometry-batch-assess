# Third-party software

This app is Apache-2.0 (see `LICENSE`). The image it builds contains the following third-party
components, installed from PyPI at the pinned versions in `requirements.txt`; their licences
apply to the image. Nothing here is GPL, so the image as a whole can be redistributed under
these permissive licences.

| Component | What | Licence | Source |
|---|---|---|---|
| FlowIO 1.3.0 | Reads FCS files (header, channels, events). | BSD-3-Clause | https://github.com/whitews/FlowIO |
| NumPy 2.2.6 | Arrays, histograms. | BSD-3-Clause | https://numpy.org |
| SciPy 1.15.3 | `scipy.stats.wasserstein_distance` (EMD). | BSD-3-Clause | https://scipy.org |
| pandas 2.2.3 | Reading samples.csv / panel.csv. | BSD-3-Clause | https://pandas.pydata.org |
| Matplotlib 3.10.3 | PNG figures. | Matplotlib licence (PSF-based, BSD-compatible) | https://matplotlib.org |
| contourpy, cycler, fonttools, kiwisolver, packaging, pillow, pyparsing, python-dateutil, pytz, six, tzdata | Dependencies of the above (versions in `requirements.txt`). | BSD-3-Clause / BSD / MIT / MIT-CMU (Pillow) / Apache-2.0 or BSD (packaging, python-dateutil) / Apache-2.0 (tzdata) | PyPI |
| awslambdaric 2.2.1 (+ simplejson 4.2.0) | AWS Lambda runtime interface client, for the Lambda runtime. | Apache-2.0 (simplejson: MIT / AFL-2.1) | https://github.com/aws/aws-lambda-python-runtime-interface-client |
| python:3.12.11-slim-bookworm | Base image (Debian 12 + CPython 3.12). | PSF (Python); Debian packages under their own free licences | https://hub.docker.com/_/python |

## Not bundled

- **CytoBatchFlagR** (https://github.com/AndorfLab/CytoBatchFlagR, GPL-3.0) does the same kind
  of assessment in R. No code was copied from it; this app's method is described in the README
  and implemented from the papers cited there.
- **Demo data** (not in the image or the repository): Koenig Z, Andorf S, Tilburgs T, Eswar S.
  CytoBatchFlagR: Spectral flow cytometry dataset. Zenodo 2025. doi:10.5281/zenodo.15388817,
  CC BY 4.0.
