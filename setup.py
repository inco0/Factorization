from distutils.core import setup
from Cython.Build import cythonize

modules_to_be_cythonized = ["fermat_square_difference.pyx", "trial_division.pyx"]


def add_package_to_modules(module: str) -> str:
    return "factorization_algorithms/" + module


setup(ext_modules=cythonize(map(add_package_to_modules, modules_to_be_cythonized),
                            language_level="3"))
