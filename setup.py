from setuptools import setup, find_packages 

VERSION = '0.0.8'
DESCRIPTION = 'Associating Gene Expression to Neuroimaging Traits Pipeline (AGENT-P)'

setup(
    name='agentp', 
    version=VERSION, 
    author='Nhung Hoang', 
    author_email='nhung.hoang@vanderbilt.edu', 
    url='https://github.com/nhunghoang/neurogen',
    description=DESCRIPTION, 
    license='MIT',
    python_requires=">=3.8",
    classifiers=[
        'Programming Language :: Python :: 3', 
        ],
    packages=find_packages(),
    install_requires=[
        'bgen_reader',
        'h5py',
        'numpy',
        'pandas',
        'patsy',
        'bitarray',
        'scipy',
        'statsmodels',
        'pyliftover', 
        'setuptools',
        ], 
    package_data={
        'agentp': ['agentp/external/MetaXcan/**/*', 
                   'agentp/external/fusion_twas/**/*',
                   'agentp/external/ldsc/**/*',
                   'agentp/external/plink2R/**/*',
                   'agentp/models/*',
                  ]
    }, 
    include_package_data=True,
) 
