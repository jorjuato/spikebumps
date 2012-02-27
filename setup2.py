from distutils.core import setup
from distutils.extension import Extension
from Pyrex.Distutils import build_ext

setup(name='prime', ext_modules=[Extension("prime", ["prime.pyx"])],\
				cmdclass = {'build_ext': build_ext}       )
