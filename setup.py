from setuptools import setup, find_namespace_packages
setup(
    name = "pie",
    package_dir={"": "pie"},
    packages=find_namespace_packages(where="pie"),
    package_data={
        "pie.module": ["*.py"],
        "pie.module.lib": ["hamming.so"],
    },
    version = "0.7.1",
    description = "Phasing all-in-one evaluator",
    author = "Jun Mencius",
    author_email = "zjmeng22@m.fudan.edu.cn",
    url = "https://github.com/JMencius/pie",
    keywords = ["pie", "haplotype", "phasing", "benchmark"],
    python_requires = ">=3.7",
    install_requires = [
        "click>=8.1.8",
        "cyvcf2>=0.31.1",
        "sortedcontainers>=2.4.0",
        "pyfastx>=2.2.0",
        ],
    extras_require = {
        "dev": ["pytest", "lmdb>=1.6.2"],
        },
    entry_points={
    "console_scripts": [
        "pie = pie.pie:main",
        ],
    },
)
