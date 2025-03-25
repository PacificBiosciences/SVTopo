import re

from setuptools import setup

with open("svtopovz/__init__.py", "r") as fd:
    version = re.search(
        r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]', fd.read(), re.MULTILINE
    ).group(1)

with open("requirements.txt", "r") as f:
    requires = f.read().splitlines()


setup(
    name="svtopovz",
    version=version,
    description="command line tool for generating visualizations of complex structural variation",
    python_requires=">=3.10",
    author="Jonathan Belyeu",
    author_email="jbelyeu@pacificbiosciences.com",
    url="",
    packages=["svtopovz"],
    package_data={"": ["LICENSE", "readme.md"]},
    include_package_data=True,
    install_requires=requires,
    license="BSD",
    zip_safe=False,
    entry_points={"console_scripts": ["svtopovz = svtopovz.__main__:main"]},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
