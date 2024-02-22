import numpy as np
import pytest
from pathlib import Path


@pytest.fixture
def resource_dir():
    return Path(__file__).parent / "resource"


@pytest.fixture
def refids_lookup(resource_dir):
    dudes_db_f = resource_dir / "dudesdb.npz"
    dudes_db = np.load(dudes_db_f, allow_pickle=True)
    return dudes_db["refids_lookup"]
