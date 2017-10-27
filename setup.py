from distutils.core import setup
setup(name='dudes',
      author='Vitor C. Piro',
	  author_email='vitorpiro@gmail.com',
	  url='https://github.com/pirovc/dudes',
	  version='0.08',
      package_dir={'dudes': ''},
      packages=['dudes'],
      scripts=['DUDes.py','DUDesDB.py']
)
