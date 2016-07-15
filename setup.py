from setuptools import setup, find_packages

with open('requirements.txt') as f:
    required = list(f.read().splitlines())

# Parse the version from the fiona module.
with open('networker/__init__.py') as f:
    for line in f:
        if line.find("__version__") >= 0:
            version = line.split("=")[1].strip()
            version = version.strip('"')
            version = version.strip("'")
            continue

setup(
    name="networker",
    version=version,
    packages=find_packages(),
    description="Python library for planning distribution networks",
    long_description=open("README.md").read(),
    package_data={'networker': ['*.json']},
    include_package_data=True,
    url='https://github.com/SEL-Columbia/networker',
    install_requires=required,
    scripts = ['scripts/run_networker.py', 'scripts/run_networkplanner.py', 'scripts/compose.py']
)
