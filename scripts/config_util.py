# -*- coding: utf-8 -*-
from __future__ import annotations

import argparse
from pathlib import Path

import yaml


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def load_config(path: str | Path | None) -> dict:
    root = repo_root()
    cfg_path = Path(path) if path else root / "configs" / "cec2.yaml"
    if not cfg_path.is_absolute():
        cand = (root / cfg_path).resolve()
        cfg_path = cand if cand.exists() else cfg_path.resolve()
    cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8")) or {}
    cfg["_root"] = root
    cfg["_config_path"] = cfg_path
    agg = Path(cfg.get("aggregates_dir", "results/paper_locked"))
    out = Path(cfg.get("output_dir", "results"))
    cfg["aggregates_dir"] = agg if agg.is_absolute() else (root / agg)
    cfg["output_dir"] = out if out.is_absolute() else (root / out)
    cfg["output_dir"].mkdir(parents=True, exist_ok=True)
    return cfg


def add_config_arg(p: argparse.ArgumentParser) -> argparse.ArgumentParser:
    p.add_argument("--config", default="configs/cec2.yaml", help="YAML config path")
    return p
