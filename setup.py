import setuptools
 
f = open("README.md", "r", encoding="utf-8")
long_description = f.read()
f.close()
 
setuptools.setup(
    name="Xponge",
    version="alpha-test",
    author="Yijie Xia",  
    author_email="yijiexia@pku.edu.cn", 
    description="A package for building molecular dynamics inputs for SPONGE",
    long_description=long_description, 
    long_description_content_type="text/markdown",
    url="https://gitee.com/gao_hyp_xyj_admin/xponge",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3 :: Only",
        "Development Status :: 3 - Alpha",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6', 
)