from setuptools import setup, find_packages

setup(
    name="reedsolomon",
    version="0.0.1",
    author="Archange",
    author_email="archange_paradise@proton.me",
    description="A Python package for Reed-Solomon error correction codes",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Archange-py/Reed-Solomon_Codes",
    project_urls={
        "Bug Tracker": "https://github.com/Archange-py/Reed-Solomon_Codes/issues",
        "Documentation": "https://github.com/Archange-py/Reed-Solomon_Codes#readme",
        "Source Code": "https://github.com/Archange-py/Reed-Solomon_Codes",
    },
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.12.4",
    install_requires=open("requirements.txt", encoding="UTF-8").read().splitlines(),
)