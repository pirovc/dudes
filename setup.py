from distutils.core import setup
from dudes import VERSION

install_requires = ["numpy", "pandas", "pyfaidx", "tqdm", "biopython", "matplotlib"]

setup(
    name="dudes",
    description="DUDes: a top-down taxonomic profiler for metagenomics",
    author="Vitor C. Piro",
    author_email="vitorpiro@gmail.com",
    url="https://github.com/pirovc/dudes",
    version=VERSION,
    packages=["dudes"],
    scripts=["DUDes.py", "DUDesDB.py", "map_peptides.py"],
    install_requires=install_requires,
)
