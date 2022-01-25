from setuptools import setup, find_packages

from os import path

import sys
if sys.version_info < (3,9):
    sys.exit('Sorry, Python < 3.9 is not supported')

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="cdftpy",
    version="1.0.0",
    author='Marat Valiev and Gennady Chuev',
    author_email='marat.valiev@gmail.com',
    description='Classical density functional theory code',
    long_description=long_description,
    long_description_content_type='text/markdown',
    license="GPL",
    classifiers=[
              'License :: OSI Approved :: GNU General Public License (GPL)',
              'Development Status :: 4 - Beta',
              'Environment :: Console',
              'Operating System :: MacOS :: MacOS X',
              'Operating System :: Microsoft :: Windows',
              'Operating System :: POSIX',
              'Programming Language :: Python',
              ],
    packages=find_packages(exclude=('tests*',)),
    package_data={ '':['data/*','examples/cdft1d/*']},
    install_requires=[
            'scipy>=1.7.3',
            'numpy>=1.21.5',
            'click>=8.0.3',
            'holoviews>=1.14.7',
            'panel>=0.12.6',
            'prettytable>=2.5.0'
    ],
    entry_points={
        'console_scripts': [
            'rism1d = cdftpy.cdft1d.cli:rism1d_run_input',
            'rsdft1d = cdftpy.cdft1d.cli:rsdft1d_run_input',
            'cdft1d = cdftpy.cdft1d.cli:cdft_cli'
        ]
    }
)
