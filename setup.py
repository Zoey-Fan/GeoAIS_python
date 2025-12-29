###python setup.py build_ext --inplace

from setuptools import setup, Extension
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        "my_rtree",  # 模块名称
        ["query3dRtree.cpp"],  # C++ 文件路径
        include_dirs=["."],  # 包含头文件的路径
        language="c++",
    ),
]

setup(
    name="my_rtree",
    version="0.1",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
