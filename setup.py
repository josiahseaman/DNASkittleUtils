from setuptools import setup, find_packages

setup(
    name='DNASkittleUtils',
    version='1.0.6',
    description='Bioinformatics functions that have been useful in multiple projects.  Manipulating FASTA files, executing pipelines, etc.',
    author='Josiah Seaman',
    author_email='josiah.seaman@gmail.com',
    license='BSD',
    packages=find_packages(exclude=('tests', 'example')),
    include_package_data=True,
    install_requires=[],
    url='https://github.com/josiahseaman/DNASkittleUtils',
    download_url='https://github.com/josiahseaman/DNASkittleUtils',  # TODO: post a tarball
    keywords=['bioinformatics', 'dna', 'fasta'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
)