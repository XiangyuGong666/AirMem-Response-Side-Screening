"""Frozen entry point: keep CLI output safe under the Windows GUI bootloader."""
import os
import sys

from airmem_response_screening.__main__ import main

if __name__ == "__main__":
    if sys.stdout is None:
        sys.stdout = open(os.devnull, "w", encoding="utf-8")
    if sys.stderr is None:
        sys.stderr = open(os.devnull, "w", encoding="utf-8")
    raise SystemExit(main())
