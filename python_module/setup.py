# setup.py - Package configuration
from setuptools import setup, find_packages

setup(
    name="siril",
    version="0.1.32",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.20.0",
        "packaging>=21.0",
        "pywin32>=300; platform_system=='Windows'"
    ],
    author="Team free-astro",
    license="GPLv3+",
    description="Python interface for Siril astronomical image processing",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://siril.org",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Copyright :: (c) Team free-astro 2024",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
    ],
    python_requires=">=3.9",
)
