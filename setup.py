import os
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
# def read(fname):
#     return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name="complex_sim",
    version="0.0.0",
    author="Heng Ma",
    author_email="heng.ma@anl.gov",
    description=(
        """A python package that automates building topologies for 
        complex protein systems with both Amber tleap and GMX pdb2gmx. """),
    # license = "BSD",
    # keywords = "example documentation tutorial",
    # url = "http://packages.python.org/an_example_pypi_project",
    # packages=['an_example_pypi_project', 'tests'],
    # long_description=read('README'),
    # classifiers=[
    #     "Development Status :: 3 - Alpha",
    #     "Topic :: Utilities",
    #     "License :: OSI Approved :: BSD License",
    # ],
)
