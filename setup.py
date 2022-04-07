from setuptools import setup
from dsautils.version import get_git_version

setup(name='dsa110-voltage-imaging',
      version=get_git_version(),
      url='http://github.com/dsa110/dsa110-voltage-imaging/',
      author='Dana Simard',
      author_email='dana.simard@astro.caltech.edu',
      packages=['dsavim'],
      package_data={'dsavim':['data/*']},
      install_requires=['astropy',
                        'casatools',
                        'casatasks',
                        'casadata',
                        'matplotlib',
                        'numpy==1.19.5',
                        'pytest',
                        'codecov',
                        'coverage',
                        'pyyaml',
                        'scipy',
      ],
)
