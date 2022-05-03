from setuptools import setup

setup(
    name = "spSeudoMap",
    version = "0.0.1",
    description = "Cell type mapping of spatial transcriptomics using unmatched single-cell RNA-seq data",
    url = "https://github.com/bsungwoo/spSeudoMap.git",
    author = "Sungwoo Bae, Hongyoon Choi",
    install_requires = ["scanpy==1.5.1","pandas==1.3.5","numpy==1.21.6",
                        "h5py==2.10.0", "jupyter",
                        "keras==2.3.1", "tensorflow==1.14.0", "tensorflow-gpu==1.14.0"]
)