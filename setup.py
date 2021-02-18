"""A setuptools based setup module for pymagnet

"""

from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / 'README.md').read_text(encoding='utf-8')


setup(
    name='pymagnet',
    version='0.1.0',
    description='User Friendly Magnetic Field Calculations', 
    long_description=long_description,  
    long_description_content_type='text/markdown',
    url='https://github.com/pdunne/pymagnet',
    author='Peter Dunne', 
    author_email='peter.dunne@applied-magnetism.com',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3 :: Only',
    ],

    keywords='magnetism, magnetic fields',
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    python_requires='>=3.6, <4',
    install_requires=['numpy', 'matplotlib'],  # Optional


    # If there are data files included in your packages that need to be
    # installed, specify them here.
    # FIXME:
    # package_data={  # Optional
    #     'sample': ['package_data.dat'],
    # },
    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/distutils/setupscript.html#installing-additional-files
    #
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    # data_files=[('my_data', ['data/data_file'])],  # Optional



    # FIXME:
    project_urls={  # Optional
        'Bug Reports': 'https://github.com/pdunne/pymagnet/issues',
        # 'Funding': 'https://donate.pypi.org',
        # 'Say Thanks!': 'http://saythanks.io/to/example',
        'Source': 'https://github.com/pdunne/pymagnet/',
    },
)
