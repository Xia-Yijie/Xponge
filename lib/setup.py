from setuptools import setup, Extension

setup(
    name="XpongeLib",
    version="0.0.8.8",
    author="Yijie Xia",  
    author_email="yijiexia@pku.edu.cn", 
    description="the C++ Lib for the package for building molecular dynamics inputs for SPONGE",
    url="https://gitee.com/gao_hyp_xyj_admin/xponge",
    package_data = {"parmchk_mod": ["*.c", "*.h"]},
    ext_modules = [Extension("XpongeLib", ["main.c"])]
)
