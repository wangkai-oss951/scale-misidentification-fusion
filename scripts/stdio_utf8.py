# -*- coding: utf-8 -*-
"""Set UTF-8 stdout on Windows consoles."""
from __future__ import annotations

import sys


def configure_stdio() -> None:
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
