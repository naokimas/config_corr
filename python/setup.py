from setuptools import setup, find_packages
from codecs import open

setup(
    name='configcorr',
    version='0.0.1',
    description='Python code for configuration model for correlation/covariance matrices',
    long_description='Python code for configuration model for correlation/covariance matrices',
    url='https://github.com/naokimas/config_corr',
    author='Naoki Masuda',
    author_email='naoki.masuda@bristol.ac.uk',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Software Development',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    keywords='correlation matrix configuration model network',
    packages=find_packages(),
    install_requires=['numpy', 'cvxpy'],
)
