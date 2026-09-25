import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SPEC = REPO_ROOT / "data" / "site_sample.json"


def sample_from_spec(spec_path):
    spec_path = Path(spec_path).resolve()
    sample = next(iter(json.loads(spec_path.read_text())["samples"].values()))
    root_file = Path(sample["files"][0])
    if not root_file.is_file():
        raise FileNotFoundError(f"Sample file is unavailable: {root_file}. Pass --spec with a local sample spec.")
    return str(root_file), sample["trees"][0]
