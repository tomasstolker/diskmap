#!/usr/bin/env python

from setuptools import setup

try:
    from pip._internal.req import parse_requirements
except ImportError:
    from pip.req import parse_requirements

reqs = parse_requirements('requirements.txt', session='hack')
reqs = [str(ir.req) for ir in reqs]

setup(
    name='diskmap',
    version='0.0.9',
    description='Scattered light mapping of protoplanetary disks',
    long_description=open('README.rst').read(),
    long_description_content_type='text/x-rst',
    author='Tomas Stolker',
    author_email='tomas.stolker@phys.ethz.ch',
    url='https://github.com/tomasstolker/diskmap',
    packages=['diskmap'],
    package_dir={'diskmap':'diskmap'},
    include_package_data=True,
    install_requires=reqs,
    license='MIT',
    zip_safe=False,
    keywords='diskmap',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
)
