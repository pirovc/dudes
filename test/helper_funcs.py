import hashlib
import os
from pathlib import Path

RESSOURCE_DIR = Path(os.path.dirname(__file__)) / "ressource"
DUDES_DIR = Path(os.path.dirname(__file__)).parent
SAMPLEDATA_DIR = DUDES_DIR / "sampledata"


def md5sum(filename):
    data = open(filename, 'rb').read()
    return hashlib.md5(data).hexdigest()
