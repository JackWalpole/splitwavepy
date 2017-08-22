from setuptools import setup

setup(
    name = 'splitting',
    version = '',
    description = "Shearwave splitting measurement tools",
    author = 'Jack Walpole',
    author_email = 'j.walpole@gmail.com',
    license = 'MIT',
    url = 'https://github.com/JackWalpole/splitting',
    packages = ['splitting'],
    install_requires = [
        'matplotlib >= 1.1',
        'numpy >= 1.1',
        'scipy']
)