from __future__ import annotations

import subprocess
from pathlib import Path


def main() -> None:
    script_path = Path(__file__).with_name("run_full_pipeline.sh")
    raise SystemExit(subprocess.run(["bash", str(script_path)], check=False).returncode)


if __name__ == "__main__":
    main()
