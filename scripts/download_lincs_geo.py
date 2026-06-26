#!/usr/bin/env python3
"""Download LINCS L1000 metadata and Level 5 signatures from GEO."""

from __future__ import annotations

import argparse
import gzip
import shutil
from pathlib import Path
from urllib.request import urlretrieve

from tqdm import tqdm


RELEASES = {
    "GSE92742": {
        "series_dir": "GSE92nnn",
        "out_dir": "data/raw/lincs_gse92742",
        "files": {
            "level5_gctx_gz": "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz",
            "sig_info": "GSE92742_Broad_LINCS_sig_info.txt.gz",
            "gene_info": "GSE92742_Broad_LINCS_gene_info.txt.gz",
            "pert_info": "GSE92742_Broad_LINCS_pert_info.txt.gz",
            "cell_info": "GSE92742_Broad_LINCS_cell_info.txt.gz",
            "sha512": "GSE92742_SHA512SUMS.txt.gz",
        },
    },
    "GSE70138": {
        "series_dir": "GSE70nnn",
        "out_dir": "data/raw/lincs_gse70138",
        "files": {
            "level5_gctx_gz": "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz",
            "sig_info": "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz",
            "gene_info": "GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz",
            "pert_info": "GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz",
            "cell_info": "GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz",
            "sha512": "GSE70138_SHA512SUMS.txt.gz",
        },
    },
}


class DownloadProgress:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.bar: tqdm | None = None

    def __call__(self, block_num: int, block_size: int, total_size: int) -> None:
        if self.bar is None:
            total = total_size if total_size > 0 else None
            self.bar = tqdm(total=total, unit="B", unit_scale=True, desc=self.path.name)
        downloaded = block_num * block_size
        self.bar.update(max(0, downloaded - self.bar.n))
        if total_size > 0 and downloaded >= total_size:
            self.bar.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--release",
        choices=sorted(RELEASES),
        default="GSE92742",
        help="GSE92742 is the LINCS Phase I release containing THP1 compound signatures.",
    )
    parser.add_argument("--out-dir", default=None)
    parser.add_argument("--skip-level5", action="store_true", help="Download only metadata files")
    parser.add_argument("--decompress-level5", action=argparse.BooleanOptionalAction, default=True)
    return parser.parse_args()


def download(url: str, dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        print(f"Exists: {dest}")
        return
    tmp = dest.with_suffix(dest.suffix + ".part")
    if tmp.exists():
        tmp.unlink()
    print(f"Downloading {url}")
    urlretrieve(url, tmp, DownloadProgress(dest))
    tmp.rename(dest)


def gunzip(src: Path, dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        print(f"Exists: {dest}")
        return
    print(f"Decompressing {src} -> {dest}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    with gzip.open(src, "rb") as fin, open(tmp, "wb") as fout:
        shutil.copyfileobj(fin, fout, length=1024 * 1024 * 64)
    tmp.rename(dest)


def main() -> None:
    args = parse_args()
    release = RELEASES[args.release]
    files = release["files"]
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{release['series_dir']}/{args.release}/suppl"
    out_dir = Path(args.out_dir or release["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    wanted = files.copy()
    if args.skip_level5:
        wanted.pop("level5_gctx_gz")

    for filename in wanted.values():
        download(f"{base}/{filename}", out_dir / filename)

    if not args.skip_level5 and args.decompress_level5:
        gz = out_dir / files["level5_gctx_gz"]
        gunzip(gz, out_dir / gz.name.removesuffix(".gz"))


if __name__ == "__main__":
    main()
