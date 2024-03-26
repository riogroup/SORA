from setuptools import setup, find_packages
from os import path

# read the contents of README.rst
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='sora-astro',
    packages=find_packages(),
    package_data={'sora': ['data/*']},
    version='0.4.0dev',
    license='MIT',
    description='Stellar Occultation Library',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='SORA Team',
    author_email='rio.occteam@gmail.com',
    url='https://github.com/riogroup/SORA',
    keywords=['science', 'astronomy', 'occultation'],
    install_requires=[
        'numpy>=1.26',
        'pyerfa>=2.0',
        'astropy>=6.0',
        'astroquery>=0.4.7',
        'spiceypy>=6.0.0',
        'matplotlib>=3.8.3',
        'scipy>=1.12',
        'requests',
        'tqdm>=4.66',
        'shapely>=2.0.3',
        'cartopy>=0.22',
    ],
    python_requires=">=3.11, <4",
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research ',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.11',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
    project_urls={
        'Bug Reports': 'https://github.com/riogroup/SORA/issues',
        'Documentation': 'https://sora.readthedocs.io/',
    },
)
