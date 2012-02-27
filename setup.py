from distutils.core import setup
from distutils.extension import Extension
from Pyrex.Distutils import build_ext

setup(name='spikes', ext_modules=[Extension("spikes", ["spikes.pyx"])],\
				cmdclass = {'build_ext': build_ext}       )
