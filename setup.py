from setuptools import setup, find_packages

with open('requirements.txt') as f:
    required = list(f.read().splitlines())

setup(
    name="networker",
    version="0.0.8",
    packages=find_packages(),
    description="Python library for planning distribution networks",
    long_description=open("README.md").read(),
    package_data={'networker': ['*.json']}, 
    include_package_data=True,
    url='https://github.com/SEL-Columbia/networker',
    install_requires=required, 
    scripts = ['scripts/run_networker.py', 'scripts/run_networkplanner.py']
)
