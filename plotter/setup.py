import re

from setuptools import setup

with open("cnidaria_plotter/__init__.py", "r") as fd:
    version = re.search(
        r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]', fd.read(), re.MULTILINE
    ).group(1)

with open("requirements.txt", "r") as f:
    requires = f.read().splitlines()


setup(
    name="cnidaria_plotter",
    version=version,
    description="command line tool for generating vizualizations of complex structural variation",
    author="Jonathan Belyeu",
    author_email="jbelyeu@pacificbiosciences.com",
    url="",
    packages=["cnidaria_plotter"],
    package_data={"": ["LICENSE", "readme.md"]},
    include_package_data=True,
    install_requires=requires,
    license="BSD",
    zip_safe=False,
    entry_points={
        "console_scripts": ["cnidaria_plotter = cnidaria_plotter.__main__:main"]
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
