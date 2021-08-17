from setuptools import setup, find_packages

setup(
    name='spaceflight-playground',
    version="0.0.1",
    setup_requires=['setuptools'],
    use_scm_version=False,
    install_requires=[
        'numpy',
        'casadi',
        'matplotlib',
    ],
    python_requires='>=3.6',
    package_dir={'': 'src'},
    packages=find_packages('src'),
    license="MIT",
    author="Paul Daum",
    author_email="paul.daum@posteo.de",
    description="Models, Optimal Control Problems and other spaceflight-related stuff",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    entry_points={
        'console_scripts': [
            #  'start = ocp1.main',
        ]
    },
)