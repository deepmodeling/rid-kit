# -*- coding: utf-8 -*-

from os import path
import setuptools
import datetime
from pathlib import Path

today = datetime.date.today().strftime("%b-%d-%Y")

# define constants
INSTALL_REQUIRES = (Path(__file__).parent / "requirements.txt").read_text().splitlines()

readme_file = Path(__file__).parent / "README.md"
readme = readme_file.read_text(encoding="utf-8")

setuptools.setup(
    name='rid',
    author="Yanze Wang, Jiahao Fan",
    use_scm_version={'write_to': 'rid/_version.py'},
    author_email="yanze039@mit.edu;jiahaofan@pku.edu.cn",
    description="RiD package for enhanced sampling",
    setup_requires=['setuptools_scm'],
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/deepmodeling/rid-kit",
    python_requires=">=3.6",
    packages=[
        "rid",
        "rid/common",
        "rid/common/gromacs",
        "rid/common/lammps",
        "rid/common/sampler",
        "rid/common/plumed",
        "rid/entrypoint",
        "rid/flow",
        "rid/nn",
        "rid/op",
        "rid/select",
        "rid/superop",
        "rid/task",
        "rid/utils",
        "rid/tools",
        "rid/common/tensorflow",
    ],
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
    ],
    package_data={'rid/template': ['*.json', '*.sh', '*.mdp']},
    include_package_data=True,
    keywords='enhanced sampling reinforced dynamics RiD',
    install_requires=INSTALL_REQUIRES,

    entry_points={
        "console_scripts": [
            "rid = rid.entrypoint.main:main"
        ]
    },
)
