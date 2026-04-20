#!/usr/bin/env python3
from pathlib import Path
import re
import sys

ROOT = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path.cwd()

"""
Run this file inside charrm to modify cetain table sizes to adapt huge bonds lists from fluctuation matching
"""


files = {
    "source/ltm/dimens_ltm.F90": [
        (
            re.compile(
                r"integer,parameter\s*::\s*MAXATC\s*=\s*1400,\s*MAXCB\s*=\s*3000,\s*MAXCH\s*=\s*6400,\s*MAXCI\s*=\s*1200,\s*&\s*\n\s*MAXCP\s*=\s*20000,\s*MAXCT\s*=\s*50000\s*!,\s*MAXITC\s*=\s*500",
                re.MULTILINE,
            ),
            "integer,parameter :: MAXATC = 10000, MAXCB = 500000, MAXCH = 6400, MAXCI = 1200, &\n"
            "         MAXCP = 20000, MAXCT = 50000 !, MAXITC = 500",
        ),
        (
            re.compile(
                r"#if KEY_BLOCK==1 /\*ldm\*/\s*\n\s*integer,parameter\s*::\s*IATBMX\s*=\s*32\s*! RLH - 16 wasn't enough\s*\n#else /\*\*/\s*\n\s*integer,parameter\s*::\s*IATBMX\s*=\s*8\s*\n#endif",
                re.MULTILINE,
            ),
            "#if KEY_BLOCK==1 /*ldm*/\n"
            "    integer,parameter :: IATBMX = 1000 ! RLH - 16 wasn't enough\n"
            "#else /**/\n"
            "    integer,parameter :: IATBMX = 1000\n"
            "#endif",
        ),
    ],
    "source/io/rtfio.F90": [
        (
            re.compile(AT
                r"integer,\s*parameter\s*::\s*mxtabl\s*=\s*50",
                re.MULTILINE,
            ),
            "integer, parameter :: mxtabl=500000",
        ),
    ],
}

changed_any = False

for relpath, patches in files.items():
    path = ROOT / relpath
    if not path.exists():
        print(f"[ERROR] Missing file: {path}")
        continue

    text = path.read_text()
    original = text

    for pattern, replacement in patches:
        new_text, n = pattern.subn(replacement, text, count=1)
        if n == 0:
            print(f"[WARN] Pattern not found in {relpath}")
        else:
            text = new_text
            print(f"[OK] Patched {relpath}")

    if text != original:
        path.write_text(text)
        changed_any = True

if changed_any:
    print("\nDone. Rebuild CHARMM from a clean tree:")
    print("  ./install.com gnu clean")
    print("  rm -rf build/ exec/ lib/")
    print("  ./install.com gnu LITE")
else:
    print("\nNo changes were made.")