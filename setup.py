from setuptools import setup, find_packages

setup(
    name='mach2',
    version='0.1',
    packages=find_packages(include=['mach2', 'mach2.*']),
    install_requires=[
        'graphviz',
        'gurobi>=9.0.0,<10.0.0 @ https://conda.anaconda.org/gurobi',
        'numpy==1.24.3',
        'pandas==2.0.1',
        'pyomo==6.5.0',
        'networkx==3.1'
    ],
    entry_points={
        'console_scripts': ['mach2=mach2.mach:main']
    }
)