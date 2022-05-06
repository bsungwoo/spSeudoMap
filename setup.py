from setuptools import setup, find_packages

setup(
    name = "spSeudoMap",
    version = "1.0.0",
    description = "Cell type mapping of spatial transcriptomics using unmatched single-cell RNA-seq data",
    url = "https://github.com/bsungwoo/spSeudoMap.git",
<<<<<<< HEAD
    packages = find_packages(include=['spSeudoMap', 'spSeudoMap.*']),
=======
    packages=find_packages(include=['spSeudoMap', 'spSeudoMap.*']),
>>>>>>> f20745a9afebe7fc9a94ae592e8d9d462b9dac76
    author = "Sungwoo Bae, Hongyoon Choi",
    install_requires = ["scanpy==1.5.1","pandas==1.3.5","numpy==1.21.6",
                        "h5py==2.10.0", "jupyter",
                        "keras==2.3.1", "tensorflow==1.14.0", "tensorflow-gpu==1.14.0"]
)
