# src/deeprbp/explainability_module/differential_regulation/utils_reg.py

import re
from typing import List

def _slugify_condition(name: str) -> str:
    """
    Turn arbitrary condition labels into safe folder names.
    - Keeps alphanumerics, '_', '-', '.'
    - Collapses whitespace into '_'
    """
    name = name.strip()
    name = re.sub(r"\s+", "_", name)
    name = re.sub(r"[^A-Za-z0-9_\-\.]+", "_", name)
    name = re.sub(r"_+", "_", name).strip("_")
    return name or "Condition"

def _parse_conditions(raw: str) -> List[str]:
    """
    Parse conditions from:
      - "A,B,C"
      - "A;B;C"
      - multiple --conditions flags are NOT supported here; keep it simple.
    """
    if not raw or not raw.strip():
        return []
    # accept comma or semicolon
    parts = re.split(r"[;,]", raw.strip())
    return [p.strip() for p in parts if p.strip()]
