from setuptools import setup, find_packages

NAME = 'dnaici'
VERSION = 0.1
DESCRIPTION = 'Installation of DNAICI Package'
LONG_DESCRIPTION = 'A package to help identify intra-chromosomal communities and differential network by integrating Hi-C data and other epigenetic information'

setup(name = NAME,
      version = VERSION,
      description = DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      packages=['dnaici',
                'dnaici.preprocess',
                'dnaici.tools',
                'dnaici.analysis',
                'dnaici.parameter'],
      entry_points={'console_scripts': ['dnaici = dnaici.dnaici:main']},
      include_package_data=True,
      install_requires=['setuptools',
                        'pandas',
                        'numpy',
                        'argparse',
                        'matplotlib',
                        'seaborn',
                        'scipy'
                        ],
      author='Zhihao Yao and Junbai Wang',
      author_email='zhhyao823917@gmail.com'
      )



