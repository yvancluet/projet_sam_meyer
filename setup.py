from setuptools import setup

setup(name='TSC',
      version='0.7',
      description='Python package for the simulation of Transcription-Supercoiling Coupling in Bacteria',
      url='https://github.com/sammeyer2017/TCDS-v2',
      author='Bilal El Houdaigui',
      author_email='bilal.elhoudaigui@gmail.com',
      license='MIT',
      packages=['TSC'],
      zip_safe=False,
      install_requires=['pandas', 'numpy', 'dnaplotlib', 'matplotlib', 'scipy'])
