#!/usr/bin/env python
"""流水线启动入口。可在任意工作目录调用。

    python scripts/run_ahsi.py preflight
    python scripts/run_ahsi.py plan
    python scripts/run_ahsi.py run --run-id demo_001
    python scripts/run_ahsi.py verify
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ahsi.cli import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
