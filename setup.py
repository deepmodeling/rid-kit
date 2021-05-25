#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from os import path
# from  dpdispatcher import NAME,SHORT_CMD
import setuptools, datetime

today = datetime.date.today().strftime("%b-%d-%Y")

install_requires=["numpy", "sklearn", "dpdispatcher", "tensorflow"]

setuptools.setup(
    name='rid',
    author="RiD Team",
    author_email="",
    description="Python RiD for enhanced sampling",
    long_description="",
    long_description_content_type="text/markdown",
    url="",
    python_requires=">=3.6",
    packages=['rid', 'rid/lib', 'rid/lib/gen', 'rid/lib/nn', 'rid/template'],
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
    ],
    package_data={'rid/template': ['*.json', '*.sh', '*.mdp']},
    include_package_data=True,
    keywords='enhanced sampling reinforced dynamics RiD',
    install_requires=install_requires,    
)