# topdown-convert

Command-line tools (via `fire`) for converting Thermo RAW → mzML (inside an Apptainer/Singularity container), filtering MS/MS by FDR, and exporting MGF + PSM tables — all from an FTP directory or a single run. Designed to be run with **[uv](https://github.com/astral-sh/uv)**.

> **Entry file:** `convert.py`  
> **CLI commands:** `single-local`, `single-ftp`, `dir-ftp`

## Prerequisites

- **uv** installed (`pipx install uv` or see uv docs).
- **Apptainer/Singularity** installed and able to run containers.
- An `msconvert` container `.sif` file.  
  - Set the path with the env var `MSCONVERT_SIF_PATH`.  
  - If not set, defaults to `~/msconvert.sif`.

```bash
export MSCONVERT_SIF_PATH="/path/to/msconvert.sif"
```

## Install & Run with uv

From the project root (where `pyproject.toml` and `convert.py` live):

```bash
# (optional) sync dependencies
uv sync

# show available subcommands
uv run convert.py --help
```

Because the CLI is powered by `fire`, you call functions as subcommands with named flags.

## Outputs

* **MGF files:** `./<output_dir>/mgf-data/*.mgf`
* **PSM tables (Parquet):** `./<output_dir>/psm-tables/*.parquet`
* **Logs:** `./<output_dir>/logs/top-down-convert-YYYYDDMM.log`

## Environmental Variables

The following environmental variables can be changed to effect the behavior of the script:

- `TDC_SE_SCORE_COL`: the search engine score column, defaults to `LogPScore`
- `TDC_MSCONVERT_SIF_PATH`: the path to the msconvert container .sif file, defaults to: `~/msconvert.sif`.

## Commands

### 1) `single-local` — Filter a *local* MS run (requires local `.msf` and matching `.mzML`)

This runs only the filtering/PSM export → produces an `.mgf` next to your `.msf`.

```bash
uv run convert.py single-local \
  --msf_path /path/to/run.msf
```

**Notes**

* Expects `/path/to/run.mzML` to already exist (same basename).
* Writes a filtered `.mgf` next to the `.msf` and returns the filtered PSMs (logged).


### 2) `single-ftp` — Download & process a *single* run from FTP

Downloads `<run>.msf` and `<run>.raw` from the given FTP directory, converts RAW → mzML via Apptainer+msconvert, filters PSMs, and writes MGF/PSM outputs.

```bash
uv run convert.py single-ftp \
  --ftp_url ftp.example.org \
  --msf_file "/remote/path/to/run.msf" \
  --output_dir ./out \
  --tmp_dir_path /scratch/tmpdir \
  --ftp_username "" \
  --ftp_password "" \
  --save_psm_table true
```

**Parameters**

* `ftp_url`: FTP host (no scheme; e.g., `ftp.pride.ebi.ac.uk`)
* `msf_file`: Remote path to the `.msf` (its sibling `.raw` is expected in the same directory)
* `output_dir`: Where to place results (`mgf-data/`, `psm-tables/`, `logs/`)
* `tmp_dir_path` (optional): Where to stage downloads (defaults to system temp)
* `ftp_username`, `ftp_password` (optional)
* `save_psm_table` (default `true`): Save filtered PSMs as Parquet

### 3) `dir-ftp` — Download & process *all* `.msf` runs in an FTP directory (multi-threaded)

Lists all `.msf` files under `ftp_path`, then processes each run in parallel (downloads its `.raw`, converts, filters, writes outputs).

```bash
uv run convert.py dir-ftp \
  --ftp_url ftp.example.org \
  --ftp_path "/remote/path/to/dir" \
  --output_dir ./out \
  --tmp_dir_path /scratch/tmpdir \
  --ftp_username "" \
  --ftp_password "" \
  --n_jobs 8
```

**Parameters**

* `ftp_url`: FTP host
* `ftp_path`: Remote directory that contains matching `.msf`/`.raw` pairs
* `output_dir`, `tmp_dir_path`, `ftp_username`, `ftp_password`: as above
* `n_jobs`: Number of parallel workers (threads). Omit to run sequentially.

## Troubleshooting

* **`msconvert` not found / container path issues**
  Ensure `MSCONVERT_SIF_PATH` points to a valid `.sif` and Apptainer can run it:

  ```bash
  apptainer exec "$MSCONVERT_SIF_PATH" wine msconvert --help
  ```

* **Permissions / bind mount**
  The RAW directory must be bind-mountable (`-B <raw_dir>:/data`). The code handles this automatically but ensure the path is readable.

* **Missing `.raw` next to `.msf` on FTP**
  The pipeline expects `<basename>.msf` and `<basename>.raw` to be siblings in the same remote folder.

## Logging

Logs go to `./<output_dir>/logs/` and stream to stdout. They include:

* FTP connections and downloads
* `msconvert` stdout (merged and streamed line-by-line)
* Filtering statistics (counts after FDR/length filtering)
* Output file locations
