import argparse
import logging
from typing import List

import duckdb
import pandas as pd

from varscore.utils.logging_config import get_logger, setup_logging

if not logging.getLogger().hasHandlers():
    setup_logging(level="INFO")

logger = get_logger(__name__)

MERGE_KEYS = ["chr", "ref", "alt"]

_REF_CANDIDATES = ["REF", "ref"]
_ALT_CANDIDATES = ["ALT", "alt"]


##################
# CORE FUNCTIONS #
##################


def init_db(parquet_files: List[str]) -> duckdb.DuckDBPyConnection:
    """Build a DuckDB connection over Hive-partitioned AlphaMissense parquet files on GCS.

    The files must live under paths of the form CHROM=chrN/data_0.parquet.
    DuckDB reads the CHROM value from the directory name and prunes partitions
    at query time, so the same connection can be reused across all variant
    batches without re-scanning unneeded chromosomes.

    Args:
        parquet_files: GCS paths to AlphaMissense parquet files,
            e.g. ['gs://bucket/alphamissense/CHROM=chr1/data_0.parquet', ...].

    Returns:
        Open DuckDB connection with an 'am_data' view ready to query.
    """
    if not parquet_files:
        raise ValueError("parquet_files must not be empty")

    con = duckdb.connect()
    con.execute("INSTALL httpfs; LOAD httpfs;")
    con.execute("CREATE SECRET gcs_secret (TYPE GCS, PROVIDER CREDENTIAL_CHAIN);")

    file_list = ", ".join(f"'{f}'" for f in parquet_files)
    logger.info(f"Registering {len(parquet_files)} AlphaMissense parquet file(s)")
    con.execute(f"""
        CREATE VIEW am_data AS
        SELECT * FROM read_parquet([{file_list}], hive_partitioning = true)
    """)

    cols = [row[0] for row in con.execute("DESCRIBE am_data").fetchall()]
    logger.info(f"View columns: {cols}")

    # CHROM is synthesised from the Hive directory name; REF/ALT are in the data.
    con._am_ref_col = _detect_col(cols, _REF_CANDIDATES, "REF")
    con._am_alt_col = _detect_col(cols, _ALT_CANDIDATES, "ALT")

    return con


def score_variants(
    variants_loc: str,
    con: duckdb.DuckDBPyConnection,
) -> pd.DataFrame:
    """Query AlphaMissense scores for variants and return the merged DataFrame.

    Args:
        variants_loc: Input variants TSV (no header; columns: chr, pos, ref,
            alt, variant_id).
        con: DuckDB connection returned by init_db.

    Returns:
        variants DataFrame merged with AlphaMissense scores.
    """
    logger.info(f"Loading variants from {variants_loc}")
    variants_df = pd.read_csv(
        variants_loc,
        sep="\t",
        header=None,
        names=["chr", "pos", "ref", "alt", "variant_id"],
    )
    logger.info(f"Loaded {len(variants_df)} variants")

    ref_col = con._am_ref_col
    alt_col = con._am_alt_col

    query_keys = variants_df[["chr", "ref", "alt"]].drop_duplicates()
    con.register("query_keys", query_keys)

    # DuckDB pushes the CHROM predicate into the Hive partition filter,
    # so only the matching chromosome files are scanned.
    query = f"""
        SELECT am.*
        FROM am_data am
        INNER JOIN query_keys qk
            ON am.CHROM = qk.chr
            AND am."{ref_col}" = qk.ref
            AND am."{alt_col}" = qk.alt
    """
    logger.info("Querying AlphaMissense data...")
    am_results = con.execute(query).df()
    con.unregister("query_keys")
    logger.info(f"Retrieved {len(am_results)} matching AlphaMissense rows")

    am_results = am_results.rename(
        columns={"CHROM": "chr", ref_col: "ref", alt_col: "alt"}
    )

    merged = variants_df.merge(am_results, on=MERGE_KEYS, how="left")
    logger.info(f"Merged result: {len(merged)} rows")
    return merged


def _detect_col(cols: list, candidates: list, label: str) -> str:
    for c in candidates:
        if c in cols:
            return c
    raise ValueError(
        f"Cannot detect {label} column. Tried {candidates}. Available: {cols}"
    )


########
# MAIN #
########


def main():
    args = _parse_args()
    con = init_db(args.parquet_files)
    merged = score_variants(args.variants_loc, con)
    con.close()
    merged.to_csv(args.out_path, sep="\t", index=False)
    logger.info(f"Saved to {args.out_path}")


def _parse_args():
    parser = argparse.ArgumentParser(
        description="Annotate variants with AlphaMissense scores from pre-scored parquet data."
    )
    parser.add_argument(
        "-p", "--parquet_files", required=True, nargs="+",
        help="GCS paths to AlphaMissense parquet files "
             "(e.g. gs://bucket/alphamissense/CHROM=chr1/data_0.parquet ...).",
    )
    parser.add_argument(
        "-v", "--variants_loc", required=True,
        help="Input variants TSV (no header; columns: chr, pos, ref, alt, variant_id).",
    )
    parser.add_argument(
        "-o", "--out_path", required=True,
        help="Output TSV path for variants merged with AlphaMissense scores.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()

"""
python -m varscore.alphamissense_scoring.scoring \
    -p gs://variant-scoring-platform/alphamissense/CHROM=chr1/data_0.parquet \
       gs://variant-scoring-platform/alphamissense/CHROM=chr2/data_0.parquet \
    -v /path/to/variants.tsv \
    -o /path/to/output_with_am_scores.tsv
"""
