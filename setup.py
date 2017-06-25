#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'Click>=6.0',
    'NetCDF4',
    'scipy',
    'numpy',
    'scikit-image',
    'matplotlib',
    'mako',
    'geojson',
    'tqdm'
    # TODO: put package requirements here
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='flowmap',
    version='0.2.11',
    description="Command line utility to transform model output into a flowmap that can be used for games or gpu-based visualizations.",
    long_description=readme + '\n\n' + history,
    author="Fedor Baart",
    author_email='fedor.baart@deltares.nl',
    url='https://github.com/SiggyF/flowmap',
    packages=find_packages(),
    package_dir={
        'flowmap': 'flowmap'
    },
    entry_points={
        'console_scripts': [
            'flowmap=flowmap.cli:cli'
        ]
    },
    scripts=[
        'bin/matroos_flowmap'
    ],
    include_package_data=True,
    install_requires=requirements,
    license="GNU General Public License v3",
    zip_safe=False,
    keywords='flowmap',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6'
    ],
    test_suite='tests',
    tests_require=test_requirements
)
