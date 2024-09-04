from distutils.core import setup

from setuptools import find_packages

with open("README.md") as fh:
    long_description = fh.read()

setup(
    name="SAM_Refiner",
    version="1.4",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    description="A program for gathering variant information from mapped reads in a SAM formated files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Devon A. Gregory",
    author_email="gregoryde@missouri.edu",
    url="https://github.com/degregory/SAM_Refiner",
    license="GPL-3.0",
    entry_points={
        "console_scripts": [
            "samrefiner = samrefiner.__main__:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
)
