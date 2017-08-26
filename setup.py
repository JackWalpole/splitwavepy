from setuptools import setup

setup(
    name = 'splitwavepy',
    version = '0.1.0',
    description = "Shearwave splitting measurement tools",
    author = 'Jack Walpole',
    author_email = 'j.walpole@gmail.com',
    license = 'MIT',
    url = 'https://github.com/JackWalpole/splitwavepy',
    packages = ['splitwavepy'],
    install_requires = [
        'matplotlib >= 1.1',
        'numpy >= 1.1',
        'scipy']
)