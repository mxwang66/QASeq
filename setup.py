import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="QASeq",
    version="0.0.1",
    author="Michael X. Wang",
    author_email="xw66@rice.edu",
    description="Analysis pipeline for QASeq (Quantitative Amplicon Sequencing)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mxwang66/QASeq",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL-3.0 license",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=['QASeq'],
    python_requires=">=3.8",
)
