from pathlib import Path
from setuptools import Extension, find_packages, setup

ROOT = Path(__file__).parent

setup(
    name="pie-phasing-eval",
    version="0.11.3",

    description="Phasing all-in-one evaluator",
    long_description=(ROOT / "README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",

    author="Jun Mencius",
    author_email="zjmeng22@m.fudan.edu.cn",
    url="https://github.com/JMencius/pie",
    keywords=["pie", "haplotype", "phasing", "benchmark"],
    license_files=["LICENSE"],

    python_requires=">=3.8",

    package_dir={"": "pie"},
    packages=find_packages(where="pie"),

    ext_modules=[
        Extension(
            "pie.module.lib.hamming",
            sources=["pie/pie/module/lib/hamming.c"],
            extra_compile_args=["-O3"],
        ),
    ],

    install_requires=[
        "click>=8.1.8",
        "cyvcf2>=0.31.3",
        "sortedcontainers>=2.4.0",
        "pyfastx>=2.2.0",
        "lmdb==1.8.1; python_version < '3.9'",
        "lmdb>=2.0.0; python_version >= '3.9'",
        "numba>=0.58.1",
    ],

    extras_require={
        "test": ["pytest"],
    },

    entry_points={
        "console_scripts": [
            "pie=pie.pie:main",
        ],
    },
)

