from setuptools import setup, find_namespace_packages
setup(
    name = "pie",
    package_dir={"": "src"},
    packages=find_namespace_packages(where="src"),
    package_data={
        "longbow.model": ["*.csv"],
        "longbow.module": ["*.py"],
    },
    version = "2.3.0",
    description = "A Python program for nanopore sequencing basecalling configuration prediction",
    author = "Jun Mencius",
    author_email = "zjmeng22@m.fudan.edu.cn",
    url = "https://github.com/JMencius/pie",
    keywords = ["longbow", "ont", "configuration"],
    python_requires = ">=3.7",
    install_requires = [
        "numpy",
        "pyfastx>=2.0.2",
        "statsmodels>=0.13.5",
        "psutil"
        ],
    extras_require = {
        "dev": ["pytest"],
        },
    entry_points={
    "console_scripts": [
        "longbow = pie.pie:main",
        ],
    },
)
