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
    scripts=["bin/archipelago-simulate.py",
    #         "bin/archipelago-generate-jobs.py",
    #         "bin/archipelago-summarize-jobs.py",
    #         "bin/archipelago-summarize-trees.py",
    #         "bin/archipelago-profile-trees.py",
            ],
    url="http://pypi.python.org/pypi/archipelago/",
    license="LICENSE.txt",
    description="A Project",
    long_description=open("README.rst").read(),
    # install_requires=[ ],
)
