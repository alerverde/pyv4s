from setuptools import setup, find_packages
import os

libv4s_path = os.path.join("pyv4s", "Cpp", "libv4s.so")
if not os.path.exists(libv4s_path):
    raise FileNotFoundError(f"Couldn't find shared library in {libv4s_path}. Please compile before installation.")

setup(
    name="pyv4s",
    version="1.0",
    description="Python interface for scientists to seamlessly perform V4S calculations using high-performance C++",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Alejandro R. Verde & Nicolás A. Loubet",
    author_email="alerverde@gmail.com",
    url="https://github.com/aleverde/pyv4s",
    license="MIT",
    packages=find_packages(),
    package_data={"pyv4s": ["../pyv4s/Cpp/libv4s.so"]},
    include_package_data=True,
    install_requires=["matplotlib", "numpy", "pandas"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
