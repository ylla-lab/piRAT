from setuptools import setup, find_packages
from setuptools.extension import Extension


setup(
    name="pirat",
    version='0.1.0',
    description="Package allowing to annotate piRNAs clusters from RNA-seq datafiles",
    url='https://github.com/shuds13/pyexample',
    author='Dominik Robak',
    author_email='dominikrobak7@gmail.com',
    packages=find_packages(),  # Add this line
    entry_points={
        'console_scripts': [
            'pirat=pirat.start:everything',
        ],
    },
    install_requires=[
    ],
    test_suite='nose.collector',
    tests_require=['nose'],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
)
