#!/usr/bin/env python3
"""Download CMap2020 compound Level 5 data and metadata."""

from __future__ import annotations

import argparse
from pathlib import Path
from urllib.request import Request, urlopen

from tqdm import tqdm


BASE_URL = "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020"
FILES = {
    "sig_info": "siginfo_beta.txt",
    "gene_info": "geneinfo_beta.txt",
    "compound_info": "compoundinfo_beta.txt",
    "level5": "level5/level5_beta_trt_cp_n720216x12328.gctx",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default="data/raw/lincs_cmap2020")
    parser.add_argument("--skip-level5", action="store_true", help="Download metadata only")
    return parser.parse_args()


def download(url: str, dest: Path) -> None:
    part = dest.with_suffix(dest.suffix + ".part")
    offset = part.stat().st_size if part.exists() else 0
    headers = {"Range": f"bytes={offset}-"} if offset else {}
    request = Request(url, headers=headers)
    with urlopen(request) as response:
        resumed = offset > 0 and response.status == 206
        if offset and not resumed:
            offset = 0
        mode = "ab" if resumed else "wb"
        remaining = int(response.headers.get("Content-Length", 0))
        total = offset + remaining if remaining else None
        with open(part, mode) as fout, tqdm(
            total=total,
            initial=offset,
            unit="B",
            unit_scale=True,
            desc=dest.name,
        ) as progress:
            while chunk := response.read(8 * 1024 * 1024):
                fout.write(chunk)
                progress.update(len(chunk))
    part.replace(dest)


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    wanted = {
        key: value
        for key, value in FILES.items()
        if not (args.skip_level5 and key == "level5")
    }
    for relative in wanted.values():
        dest = out_dir / Path(relative).name
        if dest.exists() and dest.stat().st_size > 0:
            print(f"Exists: {dest}")
            continue
        download(f"{BASE_URL}/{relative}", dest)


if __name__ == "__main__":
    main()
