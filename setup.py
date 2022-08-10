import setuptools

setuptools.setup(
    name="eddymatchup",                     # This is the name of the package
    version="0.1.0",                        # The initial release version
    author="Veronica Tamsitt",              # Full name of the author
    description="eddymatchup package to match tracked mesoscale eddy database to other ocean data",
    long_description=open('README.md').read(),      # Long description read from the the readme file
    packages=[eddymatchup,],    # List of all python modules to be installed
    url='http://pypi.python.org/pypi/eddymatchup/',
    license='LICENSE.txt', 
    python_requires='>=3.6',                # Minimum version requirement of the package
    py_modules=["quicksample"],             # Name of the python package
    package_dir={'':'eddymatchup/src'},     # Directory of the source code of the package
    install_requires=[
        "numpy >= 1.20.1",
        "xarray >= 0.18.2",
    ],
)
