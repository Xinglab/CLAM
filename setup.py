#!/usr/bin/env python

from setuptools import setup
from CLAM.config import __version__

def main():
	setup(name='CLAM',
		version=__version__,
		description='CLIP-seq Analysis of Multi-mapped reads',
		author='Zijun Zhang',
		author_email='zj.z@ucla.edu',
		url='https://github.com/Xinglab/CLAM',
		packages=['CLAM', 'CLAM.stats'],
		scripts=['bin/CLAM'],
		install_requires=[
			'scipy',
			'pysam',
			'numpy',
			'statsmodels',
			'tqdm',
			'pybedtools',
			'mpmath']
		 )
	return

if __name__ == '__main__':
	main()
