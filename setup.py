from setuptools import setup

setup(
    name='sewageQualityControl',
    version='0.1',
    packages=['lib', 'test'],
    url='https://github.com/axgraf/sewageQualityControl',
    license='',
    author='Alexander Graf',
    author_email='graf@genzentrum.lmu.de',
    description='SARS-CoV-2 Sewage quality control',
    scripts=['ssqn.py'],
    install_requires=[
        'matplotlib == 3.4.3',
        'numpy == 1.23.5',
        'pandas == 2.0.2',
        'pyarrow == 12.0.1',
        'python_dateutil == 2.8.2',
        'scikit_learn == 1.2.0',
        'scipy == 1.9.3',
        'seaborn == 0.12.2',
        'tqdm == 4.65.0',
        'openpyxl == 3.0.10',
        'adjustText == 0.7.3',
        'requests == 2.25.1',
        'jwt == 2.4.0'
    ]
    #include_package_data=True,
    #package_data={'': ['data/*.dat']},
)
