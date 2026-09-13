# Figure generation

The Fortran tools write plain text and JSON rather than figures. This directory holds the plotting scripts, which are copies of those in the Python implementation, included so that the Fortran tree is self-contained. The modules they import, `dispersion.py`, `interface_mode.py`, `nsmfp.py` and `pub_style.py`, are present alongside them, so no additional path configuration is required.

```bash
python3 plotting/make_report_figures.py --datadir validation_data \
    --archive validation_data/validation_kr0_0.7.json --outdir figures

python3 plotting/plot_results.py --json validation_data/validation_kr0_0.7.json \
    --kr0 0.7 --outdir figures
```

Figures follow the house style defined in `pub_style.py`: Computer Modern typography, inward ticks on all four sides, no axis padding, and PNG, PDF and SVG output at 600 dpi.

Reading the plain-text records written by `bin/dump_amplitude` requires only `numpy`, `scipy` and `matplotlib`. Reading plotfiles directly, which `plot_results.py` does when given `--rundir`, additionally requires `yt`.

Where a Python installation is unavailable on the machine holding the output, the record files are small and may be copied elsewhere for plotting; the analysis itself has already been performed by the Fortran tools.
