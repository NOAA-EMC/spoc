#!/usr/bin/env python3
"""Minimize BUFR files by limiting subsets per subset type or category.

This script reads an input BUFR file via the ``ncepbufr`` interface and
writes a reduced BUFR file.  By default it keeps only a configurable
number of subsets for each message type found in the input.  When a
``category_mnemonic`` is provided, the limit is applied separately to
each unique value of that mnemonic.
"""

import argparse
from collections import defaultdict
from typing import Any, Tuple

import ncepbufr


def minimize_bufr(
    in_path: str,
    out_path: str,
    max_subsets: int = 1000,
    category_mnemonic: str = "",
) -> None:
    """Create a reduced BUFR file keeping only a limited number of subsets.

    Parameters
    ----------
    in_path : str
        Path to the input BUFR file.
    out_path : str
        Path where the minimized BUFR file will be written.
    max_subsets : int, optional
        Maximum number of subsets to retain for each subset type or category.
        ``1000`` by default.
    category_mnemonic : str, optional
        BUFR mnemonic used to further categorize subsets. When set, the
        ``max_subsets`` limit is applied to each unique pair of message type and
        value of this mnemonic. An empty string disables this behaviour.
    """
    subset_counts: defaultdict[Any, int] = defaultdict(int)

    bufr_in = ncepbufr.open(in_path)
    bufr_out = ncepbufr.open(out_path, "w", table=bufr_in)

    try:
        while bufr_in.advance() == 0:
            msg_type = bufr_in.msg_type
            message_open = False

            while bufr_in.load_subset() == 0:
                if category_mnemonic:
                    try:
                        cat_val = bufr_in.read_subset(category_mnemonic).squeeze()
                        if hasattr(cat_val, "filled"):
                            cat_val = cat_val.filled("")
                        if hasattr(cat_val, "flatten"):
                            cat_val = cat_val.flatten()[0]
                        key: Tuple[Any, Any] = (msg_type, str(cat_val))
                    except Exception:
                        key = (msg_type, "")
                else:
                    key = msg_type

                if subset_counts[key] >= max_subsets:
                    continue

                if not message_open:
                    if bufr_in.receipt_time and bufr_in.receipt_time > 0:
                        bufr_out.open_message(
                            msg_type,
                            bufr_in.msg_date,
                            bufr_in.receipt_time,
                        )
                    else:
                        bufr_out.open_message(msg_type, bufr_in.msg_date)
                    message_open = True

                bufr_out.copy_subset(bufr_in)
                subset_counts[key] += 1

            if message_open:
                bufr_out.close_message()
    finally:
        bufr_in.close()
        bufr_out.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reduce the size of a BUFR file by limiting subsets")
    parser.add_argument("input", help="Path to the input BUFR file")
    parser.add_argument("output", help="Path to the minimized BUFR file")
    parser.add_argument(
        "--max-subsets",
        type=int,
        default=1000,
        help="Maximum number of subsets per subset type or category (default: 1000)",
    )
    parser.add_argument(
        "--category-mnemonic",
        default="",
        help=(
            "Mnemonic used to categorize subsets. When set, the subset limit is "
            "applied separately for each message type and mnemonic value."
        ),
    )

    args = parser.parse_args()

    minimize_bufr(args.input, args.output, args.max_subsets, args.category_mnemonic)
