"""A setuptools based setup module for pymagnet
"""

from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / "pypi_README.md").read_text(encoding="utf-8")


setup(
    name="pymagnet",
    version="0.3.4",
    description="User Friendly Magnetic Field Calculations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pdunne/pymagnet",
    author="Peter Dunne",
    author_email="peter.dunne@applied-magnetism.com",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3 :: Only",
    ],
    keywords="magnetism, magnetic fields",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.6, <4",
    install_requires=["numpy", "numpy-stl", "matplotlib", "numba", "plotly"],
    project_urls={
        "Documentation": "https://pdunne.github.io/pymagnet/",
        "Bug Reports": "https://github.com/pdunne/pymagnet/issues",
        "Source": "https://github.com/pdunne/pymagnet/",
    },
)
