# src/deeprbp/util/logger.py

import os
import torch

class Logger:
    def __init__(self, verbose=1) -> None:
        if verbose < 0:
            raise ValueError("Verbose level must be a positive integer.")
        self.verbose = verbose
    def _get_rank_prefix(self) -> str:
        local_rank = os.environ.get("LOCAL_RANK", "0")
        try:
            rank_int = int(local_rank)
        except ValueError:
            return ""
        # If rank > 0 → always show
        if rank_int > 0:
            return f"[rank {rank_int}] "
        # If rank == 0 AND we have GPU → show
        if rank_int == 0 and torch.cuda.is_available():
            return "[rank 0] "
        # Otherwise (rank 0 on CPU) → no prefix
        return ""
    def log(self, msg, level=1):
        """Log a general message."""
        if level <= self.verbose:
            print(f"{self._get_rank_prefix()}{msg}", flush=True)
    def warn(self, msg, level=1):
        """Log a warning message."""
        if level <= self.verbose:
            print(f"{self._get_rank_prefix()}WARNING: {msg}", flush=True)
    def error(self, msg, exception_type=RuntimeError, level=1):
        """Log an error message and raise a specified exception."""
        if level <= self.verbose:
            print(f"{self._get_rank_prefix()}ERROR: {msg}", flush=True)
            raise exception_type(msg)