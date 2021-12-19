import setuptools
 
f = open("README.md", "r", encoding="utf-8")
long_description = f.read()
f.close()

#for formal
setuptools.setup(
    name="Xponge",
    version="0.0.8.4.2",
    author="Yijie Xia",  
    author_email="yijiexia@pku.edu.cn", 
    description="A package for building molecular dynamics inputs for SPONGE",
    long_description=long_description, 
    long_description_content_type="text/markdown",
    url="https://gitee.com/gao_hyp_xyj_admin/xponge",
    packages=setuptools.find_packages(),
    package_data = {"":['*.mol2', '*.frcmod', '*.dat', '*.itp']},
    install_requires = ["numpy", "pubchempy"],
    classifiers=[
        "Programming Language :: Python :: 3 :: Only",
        "Development Status :: 4 - Beta",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6', 
)
