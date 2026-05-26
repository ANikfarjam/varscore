import argparse
import logging
import os

import duckdb
import pandas as pd

from varscore.utils.logging_config import get_logger, setup_logging

if not logging.getLogger().hasHandlers():
    setup_logging(level="INFO")

logger = get_logger(__name__)

MERGE_KEYS = ["chr", "ref", "alt"]

# Candidate names for each column in the AlphaMissense parquet
_CHR_CANDIDATES = ["#CHROM", "CHROM", "chrom", "chr", "chromosome"]
_REF_CANDIDATES = ["REF", "ref"]
_ALT_CANDIDATES = ["ALT", "alt"]


##################
# CORE FUNCTIONS #
##################


def score_variants(
    parquet_path: str,
    variants_loc: str,
    out_path: str,
) -> None:
    """Query AlphaMissense pre-scored data for a batch of variants.

    Initialises a DuckDB in-memory instance from parquet file(s), filters rows
    matching the input variants on (chr, ref, alt), merges the AlphaMissense
    columns back onto the input variants dataframe, and writes the result as TSV.

    Args:
        parquet_path: Path to a single parquet file, a glob pattern, or a
            directory containing *.parquet files.
        variants_loc: Path to input variants TSV (columns: chr, pos, ref, alt,
            variant_id — no header).
        out_path: Path to write the merged output TSV.
    """
    logger.info(f"Loading variants from {variants_loc}")
    variants_df = pd.read_csv(
        variants_loc,
        sep="\t",
        header=None,
        names=["chr", "pos", "ref", "alt", "variant_id"],
    )
    logger.info(f"Loaded {len(variants_df)} variants")

    con = duckdb.connect()
    parquet_glob = _resolve_parquet_glob(parquet_path)
    logger.info(f"Registering AlphaMissense parquet: {parquet_glob}")
    con.execute(f"CREATE VIEW am_data AS SELECT * FROM read_parquet('{parquet_glob}')")

    cols = [row[0] for row in con.execute("DESCRIBE am_data").fetchall()]
    logger.info(f"Parquet columns: {cols}")

    chr_col = _detect_col(cols, _CHR_CANDIDATES, "chromosome")
    ref_col = _detect_col(cols, _REF_CANDIDATES, "REF")
    alt_col = _detect_col(cols, _ALT_CANDIDATES, "ALT")

    query_keys = variants_df[["chr", "ref", "alt"]].drop_duplicates()
    con.register("query_keys", query_keys)

    query = f"""
        SELECT am.*
        FROM am_data am
        INNER JOIN query_keys qk
            ON am."{chr_col}" = qk.chr
            AND am."{ref_col}" = qk.ref
            AND am."{alt_col}" = qk.alt
    """
    logger.info("Querying AlphaMissense data...")
    am_results = con.execute(query).df()
    con.close()
    logger.info(f"Retrieved {len(am_results)} matching AlphaMissense rows")

    am_results = am_results.rename(
        columns={chr_col: "chr", ref_col: "ref", alt_col: "alt"}
    )

    merged = variants_df.merge(am_results, on=MERGE_KEYS, how="left")
    logger.info(f"Merged result: {len(merged)} rows")

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    merged.to_csv(out_path, sep="\t", index=False)
    logger.info(f"Saved to {out_path}")


def _resolve_parquet_glob(parquet_path: str) -> str:
    if os.path.isdir(parquet_path):
        return os.path.join(parquet_path, "*.parquet")
    return parquet_path


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
    score_variants(args.parquet_path, args.variants_loc, args.out_path)


def _parse_args():
    parser = argparse.ArgumentParser(
        description="Annotate variants with AlphaMissense scores from pre-scored parquet data."
    )
    parser.add_argument(
        "-p", "--parquet_path", required=True,
        help="Path to AlphaMissense parquet file(s) or directory containing them.",
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
run python -m varscore.alphamissense_scoring.scoring \
    -p /path/to/AlphaMissense_hg38.parquet \
    -v /path/to/variants.tsv \
    -o /path/to/output_with_am_scores.tsv
"""
