"""Normalize line endings of the exported CSV files and check their SHA-256 hashes.

rotterdam.csv and gbsg.csv are rewritten with pandas using "\n" line endings (as in the analyses);
nwtco.csv is converted to CRLF line endings (the form written by R on Windows, used in the analyses).
The analysis scripts stop if a hash differs.
"""
import hashlib
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
EXPECTED = {
    "rotterdam.csv": "62703670d3be5d49c5476fb7131b020bdca9f1c25bb4522fd6abe8f9dab61f1b",
    "gbsg.csv": "9fa0fef0575d04d3273869be28405307f1b02993751498e8a0c844541cccd5d9",
    "nwtco.csv": "b5fa39f11d2adb8f51845dd302806da7ed09c4eecb453548a70f0692ebf77549",
}


def main():
    for name in ("rotterdam.csv", "gbsg.csv"):
        p = HERE / name
        pd.read_csv(p).to_csv(p, index=False, lineterminator="\n")
    p = HERE / "nwtco.csv"
    p.write_bytes(p.read_bytes().replace(b"\r\n", b"\n").replace(b"\n", b"\r\n"))
    ok = True
    for name, sha in EXPECTED.items():
        got = hashlib.sha256((HERE / name).read_bytes()).hexdigest()
        status = "OK" if got == sha else "MISMATCH"
        ok &= got == sha
        print(f"{status:8s} {name}  {got}")
    if not ok:
        raise SystemExit("Hash mismatch: compare the survival package version with data/README.md.")


if __name__ == "__main__":
    main()
