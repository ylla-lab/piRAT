from setuptools import setup, find_packages
from setuptools.extension import Extension


setup(
    name="pirat",
    version='0.1.1',
    description="Package allowing to annotate piRNAs clusters from RNA-seq datafiles",
    url='https://github.com/ylla-lab/piRAT',
    author='Dominik Robak',
    author_email='dominikrobak03@gmail.com',
    packages=find_packages(),  # Add this line
    entry_points={
        'console_scripts': [
            'pirat=pirat.start:everything',
        ],
    },
    install_requires=[
        'pysam>=0.22.0',
        'numpy>=1.20.0',  
        'pandas>=2.0.0',
        'kneed>=0.8.0',
        'biopython>=1.80', 
        'tqdm>=4.60.0',
        'PyPDF2>=3.0.0',
        'tabulate>=0.9.0',
        'seqlogo>=5.29.9',
        'logomaker>=0.8',
        'cairosvg>=2.7.0',
        'matplotlib>=3.5.0',
        'imgkit>=1.2.0',
        'scipy>=1.10.0',
        'matplotlib-venn>=1.1.0', 
        'seaborn>=0.13.0',
        'intervaltree>=3.1.0',
        'pyranges>=0.1.0',
        'cycler>=0.12.0'
    ],
    python_requires='>=3.8', 
    test_suite='nose.collector',
    tests_require=['nose'],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],
)
