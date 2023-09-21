#!/usr/bin/env python

import pkg_resources
import setuptools

with open('requirements.txt') as req_txt:
    parse_req = pkg_resources.parse_requirements(req_txt)
    install_requires = [str(req) for req in parse_req]

setuptools.setup(
    name='diskmap',
    version='0.2.0',
    description='Scattered light mapping of protoplanetary disks',
    long_description=open('README.rst').read(),
    long_description_content_type='text/x-rst',
    author='Tomas Stolker',
    author_email='stolker@strw.leidenuniv.nl',
    url='https://github.com/tomasstolker/diskmap',
    project_urls={'Documentation': 'https://diskmap.readthedocs.io'},
    packages=['diskmap'],
    install_requires=install_requires,
    license='MIT',
    zip_safe=False,
    keywords='diskmap',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
)
