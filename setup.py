from setuptools import setup, find_packages
from phyloward import __version__

setup(
    name='phyloward',
    version=__version__,
    description='Phyloward is a tool for phylogenetic analysis using the prokaryotic core genes.',
    url='https://github.com/hylowaker/Phyloward',
    author='JaeHeung Han',
    author_email='hylowaker@gmail.com',
    license='GNU General Public License v3.0',
    packages=find_packages(exclude=['docs', 'tests*']),
    entry_points={
        'console_scripts': [
            'phyloward = phyloward.__main__:main'
        ]
    },
    python_requires='>= 3.5',
    include_package_data=True,
    classifiers=[
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
)

