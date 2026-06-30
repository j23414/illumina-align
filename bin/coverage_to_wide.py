#!/usr/bin/env python3

import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Convert long-format reference coverage table to wide format."
    )
    parser.add_argument(
        "--top-coverage",
        required=True,
        help="Input TSV file (all-top-coverage.tsv)"
    )
    parser.add_argument(
        "--readcount",
        required=False,
        help="Input TSV file (all-readcount.tsv)"
    )
    parser.add_argument(
        "--merged",
        required=True,
        help="Output TSV file"
    )

    args = parser.parse_args()

    # Read input
    data = pd.read_csv(
        args.top_coverage,
        sep="\t",
        header=None,
        names=["sample", "segment", "reference", "numreads", "coverage", "meandepth"],
        na_filter=False,
        dtype=str
    )

    # Preserve segment order from the input file
    segment_order = data["segment"].drop_duplicates().tolist()

    # Pivot to wide format
    wide = data.pivot(
        index="sample",
        columns="segment",
        values=["reference", "coverage", "numreads", "meandepth"]
    )

    # Reorder columns using the original segment order
    wide = wide.reindex(
        columns=segment_order,
        level=1
    )

    # Flatten MultiIndex columns
    wide.columns = [
        f"{value}_{segment}"
        for value, segment in wide.columns
    ]

    # Restore sample as a column
    wide = wide.reset_index()

    # Merge readcount table if provided
    if args.readcount:
        readcounts = pd.read_csv(
            args.readcount,
            sep="\t",
            header=None,
            names=["sample", "R1_readcount", "R2_readcount"],
            dtype=str,
            na_filter=False
        )

        wide = wide.merge(
            readcounts,
            on="sample",
            how="left"
        )

        # Put the readcount columns immediately after sample
        first_cols = ["sample", "R1_readcount", "R2_readcount"]
        other_cols = [c for c in wide.columns if c not in first_cols]
        wide = wide[first_cols + other_cols]

    # Write output
    wide.to_csv(
        args.merged,
        sep="\t",
        index=False
    )


if __name__ == "__main__":
    main()