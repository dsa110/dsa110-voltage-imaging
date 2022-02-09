from setuptools import setup
from dsautils.version import get_git_version

setup(name='dsa110-T3',
      version=get_git_version(),
      url='http://github.com/dsa110/dsa110-T3/',
      author='Dana Simard',
      author_email='dana.simard@astro.caltech.edu',
      packages=['dsaT3'],
      package_data={'dsaT3':['data/*']},
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
                        'etcd3',
                        'structlog',
                        'slack',
                        'slackclient',
                        'tensorflow==2.5.0',
                        'psrqpy'
      ],
)
