from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import Cython.Compiler.Options

#directive_defaults = Cython.Compiler.Options.get_directive_defaults()
#directive_defaults['linetrace'] = True
#directive_defaults['binding'] = True

from Cython.Build import cythonize

#setup(
#    ext_modules = cythonize("cython_pie.pyx", annotate=True ,nthreads=8, extra_compile_args=["-O3"])
#)

setup(
  name = 'Test app',
  ext_modules=[
    Extension('cython_pie',
              sources=['cython_pie.pyx'],
              extra_compile_args=['-O3'],
              #define_macros=[('CYTHON_TRACE', '1')],
              language='c++')
    ],
  cmdclass = {'build_ext': build_ext}
)
