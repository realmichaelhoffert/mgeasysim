from setuptools import setup, find_packages

setup(
    name="mgeasysim",
    version="0.2.0",
    packages=find_packages(),
    install_requires=[
        # List your dependencies here
    ],
    entry_points={
        'console_scripts': [
            'mgeasysim=mgeasysim.cli:main',
        ],
    },
    author="Your Name",
    author_email="your.email@example.com",
    description="A package for easy microbial genome simulations",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/mgeasysim",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)