# Fat and Water Separation in MRI

MATLAB university bioimaging project exploring multi-echo magnitude signal models, fat/water ambiguity, multi-start parameter fitting, and PDFF and R2* mapping.

Only source code is published. Subject scans, DICOM metadata, reports, and instructor feedback are excluded. Existing educational context is retained; this repository does not imply sole authorship of every lab-derived model or method.

## What the code covers

| Script | Purpose |
| --- | --- |
| `task1-3.m` | Compare single-peak and multi-peak signal models; fit noisy simulated signals with multiple starts |
| `task7.m` | Explore separate water/fat R2* estimates over a simulated amplitude grid |
| `advanced.m` | Compare squared-error and Rician likelihood objectives on simulated signals |
| `task4.m` | Optional DICOM loading and voxel-wise fitting, with input supplied locally |

## Run the simulations

Use MATLAB R2024a or newer with Optimization Toolbox. Run from the repository root:

```matlab
run_simulations
```

The runner fixes the random seed, creates `fig/`, and exports simulation figures. It does not require or load subject data. Full simulation execution requires MATLAB and the indicated toolbox; source review alone is not a numerical reproduction.

The GitHub Actions workflow runs these synthetic experiments using MATLAB R2024b and Optimization Toolbox, then saves figures as run artifacts. It does not run the optional DICOM workflow. Check the workflow status for actual execution results.

## Optional DICOM workflow

`task4.m` additionally requires Image Processing Toolbox and Parallel Computing Toolbox. The loader is tailored to magnitude images identified by `M_FFE` tags and expects `EchoTime` and `SliceLocation`; it is not a generic DICOM importer. Supply an authorised, de-identified local dataset:

```matlab
setenv('MRI_DICOM_DIR', '/absolute/path/to/deidentified/dicom')
setenv('MRI_SLICE_INDEX', '10')
run('task4.m')
```

The slice defaults to 10 or the last available slice if fewer exist. Result images are local and ignored by Git. Do not commit scans, patient identifiers, screenshots of clinical records, or protected metadata.

## Limitations

These are educational research scripts. Multi-start consistency is not proof of a globally optimal or clinically valid estimate. Bounds, noise-region selection, vendor-specific loading, and smoothing affect results. The public code has no clinical validation claim; raw patient data and original assessment reports are intentionally absent.
