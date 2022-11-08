import setuptools
 
f = open("README.md", "r", encoding="utf-8")
long_description = f.read()
f.close()

#for formal
setuptools.setup(
    name="Xponge",
    version="1.2.6.11.7",
    author="Yijie Xia",  
    author_email="yijiexia@pku.edu.cn", 
    description="A Python package to perform pre- and post-processing of molecular simulations",
    long_description=long_description, 
    long_description_content_type="text/markdown",
    url="https://gitee.com/gao_hyp_xyj_admin/xponge",
    packages=setuptools.find_packages(),
    package_data = {"":['*.mol2', '*.frcmod', '*.dat', '*.itp']},
    install_requires = ["numpy"],
    entry_points = {
        "console_scripts": ["Xponge = Xponge.__main__:main", 
                            "Xponge.mdrun = Xponge.mdrun.__main__:main"] },
    classifiers=[
        "Programming Language :: Python :: 3 :: Only",
        'Development Status :: 5 - Production/Stable',
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6', 
)
