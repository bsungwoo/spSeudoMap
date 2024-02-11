from setuptools import setup, find_packages

setup(
    name = "spSeudoMap",
    version = "1.1.0",
    description = "Cell type mapping of spatial transcriptomics using unmatched single-cell RNA-seq data",
    url = "https://github.com/bsungwoo/spSeudoMap.git",
    packages = find_packages(include=['spSeudoMap', 'spSeudoMap.*']),
    author = "Sungwoo Bae, Hongyoon Choi",
    install_requires = ["tensorflow~=2.9.0","tensorflow-gpu~=2.9.0", 
                        "pandas~=1.4.0","numpy~=1.20.0",
                        "scanpy","leidenalg","igraph",
                        "jupyter","ply","pytest"]
)
