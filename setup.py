from pathlib import Path

from setuptools import find_packages, setup

setup(
    name='bmitools',
    version='0.1.9',
    description='Analysis tools for calcium-based brain-computer interface experiments',
    url='https://github.com/andrestrocyte/Open-CaBCI_tools',
    license='GPL-3.0-only',
    python_requires='>=3.9',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy',
        'tqdm',
        'scipy',
        'matplotlib',
        'pyyaml',
        'networkx',
        'scikit-learn',
        'pandas',
        'opencv-python',
        'parmap',
        'scikit-image',
        'seaborn',
    ],
    long_description=Path('README.md').read_text(encoding='utf-8'),
    long_description_content_type='text/markdown',
)
