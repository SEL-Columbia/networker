from setuptools import setup

with open('requirements.txt') as f:
    required = list(f.read().splitlines())

setup(
    name="networker",
    version="0.0.1",
    packages=["networker"],
    description="Python library for computing least cost spatial networks",
    long_description=open("README.md").read(),
    url='https://github.com/SEL-Columbia/networker',
    install_requires=required
)
