from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("ideal-genom-qc")
except PackageNotFoundError:
    __version__ = "0.0.0"