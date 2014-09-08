try:
    from ez_setup import use_setuptools
    use_setuptools()
except ImportError:
    pass

from setuptools import setup, find_packages

# from mpld3
def get_version(f):
    """Get the version info from the mpld3 package without importing it"""
    import ast
    with open(f) as init_file:
        module = ast.parse(init_file.read())

    version = (ast.literal_eval(node.value) for node in ast.walk(module)
               if isinstance(node, ast.Assign)
               and node.targets[0].id == "__version__")
    try:
        return next(version)
    except StopIteration:
        raise ValueError("version could not be located")


setup(name='crystal',
      version=get_version("crystal/__init__.py"),
      description="statistical models on clusters of correlated genomic data",
      packages=find_packages(),
      url="https://github.com/brentp/crystal/",
      long_description=open('README.md').read(),
      platforms='any',
      classifiers=[
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 2.7',
          ],
      keywords='bioinformatics methylation correlation',
      author='brentp',
      author_email='bpederse@gmail.com',
      license='MIT License',
      include_package_data=True,
      tests_require=['nose'],
      test_suite='nose.collector',
      zip_safe=False,
      install_requires=['numpy', 'pandas', 'aclust', 'toolshed', 'statsmodels',
          'seaborn', 'scipy', 'patsy', 'numpydoc'],
      #scripts=[],
      entry_points={
      },
  )
