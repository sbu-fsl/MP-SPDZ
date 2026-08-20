#!/usr/bin/env python3

'''
Every time we sync with upstream MP-SPDZ fork, git brings in any new files
committed to MP-SPDZ/Programs/Source. These files are often irrelevant
and clutter the directory. This script implements a program whitelist that only
keeps files listed in ALLOWED_FILES.

Please keep ALLOWED_FILES up to date if you use this script...
'''

import argparse
import subprocess
from pathlib import Path

# assumes we are running from MP-SPDZ/Scripts/
SOURCE_DIR = Path("../MP-SPDZ/Programs/Source")

# Paths are relative to SOURCE_DIR.
# Please keep in alphabetical order
# Please keep updated
ALLOWED_FILES = {
    "shared-library.mpc",
    "aes_circuit.mpc",
    "aes_ctr.py",
    "aes.mpc",
    "aes.py",
    "arithmetic.py",
    "bankers_bonus.mpc",
    "cmac.py",
    "embeddings.py",
    "gf2n.py",
    "kdf-ctr.py",
    "linalg_lists.py",
    "linalg.py",
    "lrpss.py",
    "lrss_key_gen.py",
    "lrss-old.py",
    "lrss.mpc",
    "lrss.py",
    "millionaires-server.py",
    "old_threshold_key_derivation.py",
    "old-shamir.py",
    "oram_playground.py",
    "personal_client_example.py",
    "robust_lrpss.py",
    "robust_lrss.py",
    "shamir.py",
    "shamir_key_gen.py",
    "shamir_pss.py",
    "share_io.py",
    "tkdf.py",
    "utils.py",
}

for path in SOURCE_DIR.rglob("*"):
    if path.is_file() and path.relative_to(SOURCE_DIR).as_posix() not in ALLOWED_FILES:
        print(f"Deleting {path}")
        path.unlink()
