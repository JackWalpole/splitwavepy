from setuptools import setup, find_packages

setup(
    name = 'SplitWavePy',
    version = '0.3.0',
    description = "Shearwave splitting measurement tools",
    author = 'Jack Walpole',
    author_email = 'j.walpole@gmail.com',
    license = 'MIT',
    url = 'https://github.com/JackWalpole/splitwavepy',
    packages = find_packages(),
    install_requires = ['matplotlib','numpy','scipy'],
    keywords = 'shear wave splitting seismic anisotropy seismology',
)
