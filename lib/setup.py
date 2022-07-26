from setuptools import setup, Extension, find_packages

setup(
    name="XpongeLib",
    version="1.2.5.0",
    author="Yijie Xia",  
    author_email="yijiexia@pku.edu.cn", 
    description="the C++ Lib for the package for building molecular dynamics inputs for SPONGE",
    url="https://gitee.com/gao_hyp_xyj_admin/xponge",
    packages = find_packages(),
    package_data = {"": ["parmchk_mod/*.dat", "parmchk_mod/*.DAT", "parmchk_mod/CONNECT.TPL"]},
    ext_modules = [Extension("XpongeLib.backend", ["XpongeLib/main.c"])]
)
