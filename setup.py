from distutils.core import setup

setup(name='dudes',
      description='DUDes: a top-down taxonomic profiler for metagenomics',
      author='Vitor C. Piro',
      author_email='vitorpiro@gmail.com',
      url='https://github.com/pirovc/dudes',
      version='0.08',
      packages=['dudes'],
      scripts=['DUDes.py', 'DUDesDB.py'],
      install_requires=[
          "numpy",
          "pandas"
      ]
      )
