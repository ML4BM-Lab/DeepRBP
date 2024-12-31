# logger.py

class Logger:
    def __init__(self, verbose=1) -> None:
        if verbose < 0:
            raise ValueError("Verbose level must be a positive integer.")
        else:
            self.verbose = verbose

    def log(self, msg, level=1):
        if level <= self.verbose:
            print(msg, flush=True)

    def warn(self, msg, level=1):
        if level <= self.verbose:
            print(f"WARNING: {msg}", flush=True)

    def error(self, msg, exception_type=RuntimeError, level=1):
        """Log an error message and raise a specified exception."""
        if level <= self.verbose:
            print(f"ERROR: {msg}", flush=True)
            raise exception_type(msg)