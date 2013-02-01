from setuptools import setup, find_packages
from setuptools.extension import Extension
import sys, os

version = '0.1'

c_ext = Extension("facs", define_macros = [('DEBUG', '1'), ('FIFO', '1'), ('FILE_OFFSET_BITS', '64'), ('LARGE_FILE', '1')],
                           sources = ["facs.c", "tool.c", "bloom.c", "good_build.c",
                                      "suggestions.c", "lookup8.c", "file_dir.c",
                                      "simple_check_1_ge.c", "big_query.c", "simple_remove.c"],
                           extra_compile_args = ['-fopenmp'],
                           extra_link_args=['-lgomp', '-lz'])

setup(name='facs',
      version=version,
      description="FACS bloom filter implementation",
      long_description="""FACS you""",
      ext_modules=[c_ext],
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='bloom filter probabilistic',
      author='Enze Liu, Lars Arvestad, Henrik Stranneheim, Roman Valls Guimera',
      author_email='roman@scilifelab.se',
      url='http://facs.scilifelab.se/',
      license='GPLv3',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
