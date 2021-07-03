from setuptools import setup

setup(
    name="dpdata_random",
    packages=['dpdata_random'],
    entry_points={
        'dpdata.plugins': [
            'random=dpdata_random:RandomFormat'
        ]
    },
    install_requires=['dpdata', 'numpy'],
)
