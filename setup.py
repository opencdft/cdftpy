from setuptools import setup, find_packages

from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="cdftpy",
    version="0.1.17",
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
    package_data={ "":["data/*"]},
    install_requires=[
            'scipy',
            'numpy',
            'matplotlib',
            'click',
            'prompt_toolkit',
            'holoviews',
            'panel',
            'prettytable'
    ],
    entry_points={
        'console_scripts': [
            'rism1d = cdftpy.cdft1d.cli:rism1d_run_input',
            'rsdft1d = cdftpy.cdft1d.cli:rsdft1d_run_input',
            'cdft1d = cdftpy.cdft1d.cli:cdft_cli'
        ]
    }
)
