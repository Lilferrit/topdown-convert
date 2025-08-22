import dataclasses
import datetime
import ftplib
import logging
import os
import pathlib
import re
import shutil
import sqlite3
import subprocess
import tempfile
from typing import Optional

import crema
import fire
import joblib
import pandas as pd
import pyteomics.mgf
import pyteomics.mzml
import tqdm

PathLike = pathlib.Path | str

SCORE_COLUMN = os.getenv("TDC_SE_SCORE_COL", default="LogPScore")
MSCONVERT_SIF_PATH = os.getenv("TDC_MSCONVERT_SIF_PATH", default="~/msconvert.sif")

@dataclasses.dataclass
class PsmInfo:
    precursor_charge: int
    precursor_mz: float
    sequence: str


def setup_logging(
    log_dir: Optional[PathLike] = None,
    to_file: bool = True,
    filename: Optional[str] = None,
) -> None:
    """Configure logging for the topdown-convert pipeline.

    If logging has already been configured, this function does nothing.

    Parameters
    ----------
    log_dir : PathLike, optional
        Directory where log files will be written. If None, uses the
        current working directory. Ignored if `to_file` is False.
    to_file : bool, default=True
        If True, logs are written both to a file and to the console.
        If False, logs are only streamed to the console.
    filename : str, optional
        Custom log filename. If None, defaults to
        ``top-down-convert-YYYYDDMM.log``.
    """
    root_logger = logging.getLogger()
    if root_logger.handlers:  # already configured → no-op
        return

    if filename is None:
        dt_str = datetime.datetime.now().strftime("%Y%d%m")
        filename = f"top-down-convert-{dt_str}.log"

    handlers = [logging.StreamHandler()]

    if to_file:
        if log_dir is None:
            log_dir = pathlib.Path.cwd()
        else:
            log_dir = pathlib.Path(log_dir)

        log_dir.mkdir(parents=True, exist_ok=True)
        log_path = log_dir / filename
        handlers.insert(0, logging.FileHandler(log_path, mode="w"))

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] [%(funcName)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=handlers,
    )


def get_psms(msf_file_path: PathLike) -> pd.DataFrame:
    """
    Load PSM (proteoform spectrum match) data from a Thermo MSF database.

    Parameters
    ----------
    msf_file_path : PathLike
        Path to the MSF SQLite database file.

    Returns
    -------
    pd.DataFrame
        Concatenated DataFrame of target and decoy PSMs with missing
        ``ParentProteinAccessions`` values filled as empty strings.
    """
    with sqlite3.connect(msf_file_path) as con:
        TargetProteoformSpectrumMatchs = pd.read_sql(
            "SELECT * FROM TargetProteoformSpectrumMatchs", con
        )
        DecoyProteoformSpectrumMatchs = pd.read_sql(
            "SELECT * FROM DecoyProteoformSpectrumMatchs", con
        )

        TargetProteoformSpectrumMatchs = TargetProteoformSpectrumMatchs.merge(
            pd.read_sql(
                """
                SELECT
                    TargetProteoformSpectrumMatchsID,
                    MSnSpectrumInfoSpectrumID
                FROM
                    TargetProteoformSpectrumMatchsMSnSpectrumInfo
                """,
                con,
            ),
            left_on="ID",
            right_on="TargetProteoformSpectrumMatchsID",
        )

        DecoyProteoformSpectrumMatchs = DecoyProteoformSpectrumMatchs.merge(
            pd.read_sql(
                """
                SELECT
                    DecoyProteoformSpectrumMatchsID,
                    MSnSpectrumInfoSpectrumID
                FROM
                    DecoyProteoformSpectrumMatchsMSnSpectrumInfo
                """,
                con,
            ),
            left_on="ID",
            right_on="DecoyProteoformSpectrumMatchsID",
        )

    psm_df = pd.concat(
        [TargetProteoformSpectrumMatchs, DecoyProteoformSpectrumMatchs], 
        ignore_index=True
    )
    
    # The decoy tables doesn't populate this column
    psm_df["is_target"] = ~psm_df["CScore"].isna()
    return psm_df


def get_target_ms_info(
    target_df: pd.DataFrame, msf_file_path: PathLike
) -> pd.DataFrame:
    """
    Join target PSMs with spectrum-level MS information.

    Parameters
    ----------
    target_df : pd.DataFrame
        Target-only PSM DataFrame.
    msf_file_path : PathLike
        Path to the MSF SQLite database file.

    Returns
    -------
    pd.DataFrame
        DataFrame of target PSMs with additional spectrum-level metadata merged.
    """
    with sqlite3.connect(msf_file_path) as con:
        ms_info_df = pd.read_sql("SELECT * FROM MSnSpectrumInfo", con)

    target_df = target_df.merge(
        ms_info_df,
        left_on="MSnSpectrumInfoSpectrumID",
        right_on="SpectrumID",
        suffixes=["_psm_info", "_spectrum_info"],
    )

    return target_df


def get_filtered_targets(
    psm_df: pd.DataFrame,
    max_fdr: float = 0.01,
    min_pep_len: Optional[int] = 5,
    max_pep_len: Optional[int] = None,
) -> pd.DataFrame:
    """
    Apply FDR and (optional) peptide length filtering to PSMs.

    Parameters
    ----------
    psm_df : pd.DataFrame
        DataFrame containing target and decoy PSMs.
    max_fdr : float, default=0.01
        Maximum acceptable FDR (q-value).
    min_pep_len : int, optional, default=5
        Minimum peptide sequence length to retain. If None, no minimum is
        applied.
    max_pep_len : int, optional, default=None
        Maximum peptide sequence length to retain. If None, no maximum i
        applied.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame of target PSMs that pass FDR and length criteria.
    """
    psm_df = psm_df.reset_index(drop=True)
    psm_df["ParentProteinAccessions"] = psm_df["ParentProteinAccessions"].fillna("")
    psms = crema.read_txt(
        txt_files=psm_df,
        target_column="is_target",
        spectrum_columns="MSnSpectrumInfoSpectrumID",
        score_columns=SCORE_COLUMN,
        peptide_column="Sequence",
        protein_column="ParentProteinAccessions",
        protein_delim=",",
    )

    confidence = psms.assign_confidence(pep_fdr_type="psm-only", threshold="q-value")
    crema_df = confidence.confidence_estimates["psms"]
    crema_df = crema_df[crema_df["crema q-value"] <= max_fdr]
    
    if min_pep_len is not None:
        crema_df = crema_df[crema_df["Sequence"].str.len() >= min_pep_len]
        
    if max_pep_len is not None:
        crema_df = crema_df[crema_df["Sequence"].str.len() <= max_pep_len]
        
    filtered_df = psm_df.loc[crema_df.index]
    return filtered_df


def run_ms_convert(raw_file: PathLike) -> pathlib.Path:
    """
    Convert a Thermo RAW file to mzML format using msconvert inside Apptainer.

    Invokes the ``msconvert`` command wrapped in a Singularity/Apptainer
    container. The RAW file is bind-mounted into the container and converted to
    a compressed, peak-picked mzML file. Assumes the container file is at the
    path ~/msconvert.sif if the envar MSCONVERT_SIF_PATH is not set.

    Parameters
    ----------
    raw_file : PathLike
        Path to the Thermo RAW file to convert.

    Returns
    -------
    pathlib.Path
        Path to the generated mzML file (same base name as the RAW file).
    """
    raw_file = pathlib.Path(raw_file)
    raw_file_dir = raw_file.parent
    mzml_file = raw_file.with_suffix(".mzML")
    sif_path = pathlib.Path(MSCONVERT_SIF_PATH).expanduser().resolve()
    logging.info(f"Converting {raw_file.name} to {mzml_file.name}")

    cmd = [
        "apptainer", "exec",
        "-B", f"{raw_file_dir}:/data",
        str(sif_path),
        "wine", "msconvert",
        f"/data/{raw_file.name}",
        "--outdir", "/data",
        "--verbose",
        "--zlib",
        "--mzML",
        "--64",
        "--filter", "peakPicking true 1-",
    ]

    logging.info("Running command: %s", " ".join(cmd))
    process = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1
    )

    # Stream each line into logging
    with process.stdout:
        for line in process.stdout:
            logging.info(line.strip())

    retcode = process.wait()
    if retcode != 0:
        raise subprocess.CalledProcessError(retcode, cmd)

    return mzml_file


def get_filtered_mgf(msf_path: PathLike) -> pd.DataFrame:
    """
    Generate an MGF file filtered by high-confidence PSMs.

    Parameters
    ----------
    msf_path : PathLike
        Path to the MSF SQLite database file.

    Returns
    -------
    pd.DataFrame
        Filtered PSM DataFrame used to generate the MGF file.
    """
    setup_logging(to_file=False)
    
    msf_path = pathlib.Path(msf_path)
    mzml_path = msf_path.with_suffix(".mzML")
    mgf_path = msf_path.with_suffix(".mgf")
    raw_path = msf_path.with_suffix(".raw")
    
    if not mzml_path.is_file():
        if not raw_path.is_file():
            raise FileNotFoundError(f"No spectrum file found for {msf_path}")
        run_ms_convert(raw_path)

    logging.info(f"Filtering MS run: {msf_path}")
    filtered_df = get_psms(msf_path)
    filtered_df = get_filtered_targets(filtered_df)
    filtered_df = get_target_ms_info(filtered_df, msf_path)
    logging.info(f"Number of remaining MS2 scans post filtering: {len(filtered_df)}")

    filtered_scans = {
        scan_num: PsmInfo(
            precursor_charge=precursor_charge, precursor_mz=precursor_mz, sequence=seq
        )
        for scan_num, precursor_mz, precursor_charge, seq in zip(
            filtered_df["FirstScan"],
            filtered_df["MassOverCharge_spectrum_info"],
            filtered_df["OriginalPrecursorCharge_spectrum_info"],
            filtered_df["Sequence"],
        )
    }

    logging.info(f"Reading mzML file: {mzml_path}")
    mgf_spectra = list()
    scan_pattern = re.compile(r"scan=(\d+)")
    with pyteomics.mzml.read(str(mzml_path)) as mzml_spectra:
        for spectrum in tqdm.tqdm(mzml_spectra, desc="Reading mzML"):
            if spectrum["ms level"] != 2:
                continue

            scan_id = spectrum["id"]
            scan_id_match = scan_pattern.search(scan_id)

            if scan_id_match:
                scan_number = int(scan_id_match.group(1))
            else:
                continue

            if scan_number not in filtered_scans:
                continue

            psm_info = filtered_scans[scan_number]
            mgf_spectra.append(
                {
                    "params": {
                        "title": spectrum["id"],
                        "seq": psm_info.sequence,
                        "charge": f"{psm_info.precursor_charge}+",
                        "pepmass": psm_info.precursor_mz,
                    },
                    "m/z array": spectrum["m/z array"],
                    "intensity array": spectrum["intensity array"],
                }
            )

    logging.info(f"Writing MGF file: {mgf_path}")
    pyteomics.mgf.write(
        tqdm.tqdm(mgf_spectra, desc="writing mgf"), output=str(mgf_path)
    )
    return filtered_df


def get_ms_run_ftp(
    ftp_url: str,
    msf_file: PathLike,
    output_dir: PathLike,
    tmp_dir_path: Optional[PathLike] = None,
    ftp_username: str = "",
    ftp_password: str = "",
    save_psm_table: bool = True,
) -> Optional[pd.DataFrame]:
    """
    Download and process a single mass spectrometry run from an FTP server.

    Parameters
    ----------
    ftp_url : str
        Hostname of the FTP server.
    msf_file : PathLike
        Path to the `.msf` file on the FTP server.
    output_dir : PathLike
        Local directory where `.mgf` files and PSM tables are written.
    tmp_dir_path : PathLike, optional
        Directory for temporary storage of downloaded files.
        Defaults to the system temp directory.
    ftp_username : str, default=""
        Username for FTP authentication.
    ftp_password : str, default=""
        Password for FTP authentication.
    save_psm_table : bool, default=True
        Whether to save the PSM info table to disk as Parquet.

    Returns
    -------
    pandas.DataFrame or None
        The extracted PSM table if conversion succeeds, otherwise None.
    """
    if tmp_dir_path is None:
        tmp_dir_path = tempfile.gettempdir()

    output_dir = pathlib.Path(output_dir)
    tmp_dir_path = pathlib.Path(tmp_dir_path)
    msf_file = pathlib.Path(msf_file)
    ftp_path = msf_file.parent
    psm_data_table = None

    mgf_dir = output_dir / "mgf-data"
    mgf_dir.mkdir(parents=True, exist_ok=True)
    logging.info(f"Writing MGF file to: {mgf_dir}")

    with ftplib.FTP(ftp_url) as ftp:
        ftp.login(user=ftp_username, passwd=ftp_password)
        logging.info("Logged into FTP server successfully.")

        ftp.cwd(str(ftp_path))
        logging.info(f"Changed directory to: {ftp_path}")

        try:
            msf_file = pathlib.Path(msf_file)
            raw_file = msf_file.with_suffix(".raw")
            logging.info(
                f"Processing MS run: {msf_file.name}, expecting RAW file: "
                f" {raw_file.name}"
            )

            folder_name = msf_file.stem
            curr_tmp_dir = tmp_dir_path / folder_name
            curr_tmp_dir.mkdir(parents=True, exist_ok=True)

            for file in [msf_file, raw_file]:
                logging.info(f"Downloading {file.name} from FTP to {curr_tmp_dir}")
                with open(curr_tmp_dir / file.name, "wb") as f:
                    ftp.retrbinary(f"RETR {file.name}", f.write)

            raw_file = curr_tmp_dir / raw_file.name
            msf_file = curr_tmp_dir / msf_file.name
            logging.info(f"Running msconvert on {raw_file}")
            run_ms_convert(raw_file)

            logging.info(f"Filtering and extracting MGF for {msf_file}")
            psm_data_table = get_filtered_mgf(msf_file)

            mgf_file = msf_file.with_suffix(".mgf")
            shutil.copy(mgf_file, mgf_dir / mgf_file.name)
            logging.info(f"Copied MGF file to {mgf_dir / mgf_file.name}")
        except Exception:
            logging.exception(f"Failed to convert MS run {msf_file}")
        finally:
            shutil.rmtree(curr_tmp_dir)

    if not save_psm_table:
        return psm_data_table

    if psm_data_table is None:
        logging.warning("Not saving PSM info table due to conversion error")
    else:
        data_table_dir = output_dir / "psm-tables"
        data_table_dir.mkdir(exist_ok=True, parents=True)
        data_table_path = data_table_dir / msf_file.with_suffix(".parquet").name
        logging.info(f"Saving PSM info table to {data_table_path}")
        psm_data_table.to_parquet(data_table_path)

    return psm_data_table


def get_ms_directory_ftp(
    ftp_url: str,
    ftp_path: PathLike,
    output_dir: PathLike,
    tmp_dir_path: Optional[PathLike] = None,
    ftp_username: str = "",
    ftp_password: str = "",
    n_jobs: Optional[int] = None,
) -> None:
    """
    Download and process all MSF-based runs from a directory on an FTP server.

    Parameters
    ----------
    ftp_url : str
        Hostname of the FTP server.
    ftp_path : PathLike
        Remote directory on the FTP server containing `.msf` and `.raw` files.
    output_dir : PathLike
        Local directory where `.mgf` files, PSM tables, and logs are written.
    tmp_dir_path : PathLike, optional
        Local temporary directory for intermediate files. Defaults to the
        system temp directory.
    ftp_username : str, default=""
        Username for FTP authentication.
    ftp_password : str, default=""
        Password for FTP authentication.
    n_jobs : int, optional
        Number of parallel jobs to use for downloading and processing files.
        If None, runs sequentially.

    Returns
    -------
    None
    """
    output_dir = pathlib.Path(output_dir)
    log_dir = output_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(log_dir)

    logging.info(f"Connecting to FTP server: {ftp_url}")
    with ftplib.FTP(ftp_url) as ftp:
        ftp.login(user=ftp_username, passwd=ftp_password)
        logging.info("Logged into FTP server successfully.")

        ftp.cwd(str(ftp_path))
        logging.info(f"Changed directory to: {ftp_path}")
        msf_files = [f for f in ftp.nlst() if f.lower().endswith(".msf")]

    logging.info(f"Found {len(msf_files)} MSF file(s) under {ftp_path}")
    for msf_file in msf_files:
        logging.info(f"  {msf_file}")

    joblib.Parallel(n_jobs=n_jobs, prefer="threads", batch_size=1)(
        joblib.delayed(get_ms_run_ftp)(
            ftp_url,
            ftp_path / msf_file,
            output_dir,
            tmp_dir_path=tmp_dir_path,
            ftp_username=ftp_username,
            ftp_password=ftp_password,
        )
        for msf_file in msf_files
    )


def main():
    """CLI Entry"""
    fire.Fire(
        {
            "single-local": get_filtered_mgf,
            "single-ftp": get_ms_run_ftp,
            "dir-ftp": get_ms_directory_ftp,
        }
    )


if __name__ == "__main__":
    main()
