# -*- coding: utf-8 -*-

import os
import subprocess
from setuptools import setup, find_packages

# Generate version.py
__version__ = None
with open('zap/version.py') as f:
    exec(f.read())

# If the version is not stable, we can add a git hash to the __version__
if '.dev' in __version__:
    # Find hash for __githash__ and dev number for __version__ (can't use hash
    # as per PEP440)
    command_hash = 'git rev-list --max-count=1 --abbrev-commit HEAD'
    command_number = 'git rev-list --count HEAD'

    try:
        commit_hash = subprocess.check_output(command_hash, shell=True)\
            .decode('ascii').strip()
        commit_number = subprocess.check_output(command_number, shell=True)\
            .decode('ascii').strip()
    except Exception:
        pass
    else:
        # We write the git hash and value so that they gets frozen if installed
        with open(os.path.join('zap', '_githash.py'), 'w') as f:
            f.write("__githash__ = \"{}\"\n".format(commit_hash))
            f.write("__dev_value__ = \"{}\"\n".format(commit_number))

        # We modify __version__ here too for commands such as egg_info
        __version__ += commit_number


setup(
    name='zap',
    version=__version__,
    description=('ZAP (the Zurich Atmosphere Purge) is a high precision sky'
                 ' subtraction tool'),
    author='Kurt Soto',
    author_email='kurt.soto@phys.ethz.ch',
    url='https://github.com/ktsoto/zap',
    license='MIT',
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    install_requires=['numpy', 'scipy', 'astropy'],
    extras_require={'plot': ['matplotlib']},
    entry_points={
        'console_scripts': ['zap = zap.__main__:main']
    },
)
