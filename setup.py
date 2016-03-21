#! /usr/bin/env python

#from distutils.core import setup
from setuptools import setup

setup(
    name="archipelago",
    version="0.1.0",
    author="Jeet Sukumaran",
    author_email="jeetsukumaran@gmail.com",
    packages=[
        "archipelago",
        # "test"
        ],
    # these directories,
    package_data = {
            'archipelago': ['R/*.R'],
            'archipelago': ['libexec/*.R'],
        },
    include_package_data = True,
    scripts=[
            "bin/archipelago-classify.py",
            "bin/archipelago-profile-trees.py",
            "bin/archipelago-simulate.py",
            "bin/archipelago-summarize.py",
            "bin/archipelago-generate-data-files-from-tip-labels.py",
            "bin/archipelago-encode-trees.py",
            ],
    url="http://pypi.python.org/pypi/archipelago/",
    license="LICENSE.txt",
    description="A Project",
    long_description=open("README.md").read(),
    # install_requires=[ ],
)
