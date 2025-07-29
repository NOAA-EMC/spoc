# BUFR Minimizer

This utility reduces the size of an existing BUFR file by limiting the number
of subsets kept for each subset type.  It can be useful when generating small
sample BUFR files for unit tests or debugging.

## Usage

```bash
python tools/bufr_minimizer/bufr_minimizer.py INPUT_FILE OUTPUT_FILE [--max-subsets N] [--category-mnemonic MNEM]
```

- `INPUT_FILE` – path to the original BUFR file.
- `OUTPUT_FILE` – location where the minimized BUFR file will be written.
- `--max-subsets` – optional maximum number of subsets to keep for each subset
  type or category. Defaults to `1000` if not provided.
- `--category-mnemonic` – optional mnemonic used to further categorize subsets.
  When set, the subset limit applies to each unique combination of message type
  and mnemonic value.

The script relies on the `ncepbufr` Python interface provided by the
`NCEPLIBS-bufr` project.

Make sure this package is installed and available on your `PYTHONPATH` before
running the script.
