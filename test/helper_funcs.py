import hashlib
import os
from contextlib import contextmanager
from pathlib import Path

RESOURCE_DIR = Path(os.path.dirname(__file__)) / "resource"
DUDES_DIR = Path(os.path.dirname(__file__)).parent
SAMPLEDATA_DIR = DUDES_DIR / "sampledata"


def md5sum(filename):
    data = open(filename, 'rb').read()
    return hashlib.md5(data).hexdigest()


@contextmanager
def set_directory(path: Path):
    """Sets the cwd within the context

    Args:
        path (Path): The path to the cwd

    Yields:
        None
    """

    origin = Path().absolute()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(origin)
