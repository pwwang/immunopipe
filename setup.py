# -*- coding: utf-8 -*-
from setuptools import setup

packages = \
['immunopipe', 'immunopipe.scripts.TCR-counts']

package_data = \
{'': ['*'],
 'immunopipe': ['reports/*',
                'scripts/*',
                'scripts/Seurat-0/*',
                'scripts/utils/*']}

install_requires = \
['pipen',
 'pipen-args',
 'pipen-report',
 'pipen-verbose',
 'pyparam',
 'rich>=10.0.0,<11.0.0']

entry_points = \
{'console_scripts': ['immunopipe = immunopipe:main']}

setup_kwargs = {
    'name': 'immunopipe',
    'version': '0.0.0',
    'description': 'Integrative analysis for scTCR- and scRNA-seq data',
    'long_description': None,
    'author': 'pwwang',
    'author_email': 'pwwang@pwwang.com',
    'maintainer': None,
    'maintainer_email': None,
    'url': None,
    'packages': packages,
    'package_data': package_data,
    'install_requires': install_requires,
    'entry_points': entry_points,
    'python_requires': '>=3.7.1,<4.0.0',
}


setup(**setup_kwargs)

# This setup.py was autogenerated using poetry.
