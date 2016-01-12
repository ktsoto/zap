from setuptools import setup, find_packages

setup(
    name='zap',
    version='0.6-dev',
    description=('ZAP (the Zurich Atmosphere Purge) is a high precision sky'
                 ' subtraction tool'),
    author='Kurt Soto',
    author_email='sotok@phys.ethz.ch',
    # url='',
    # license='',
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    install_requires=['numpy', 'scipy', 'astropy', 'joblib'],
    extras_require={'plot': ['matplotlib']},
)
