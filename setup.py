from setuptools import setup

setup(
    name='sewageQualityControl',
    version='0.1',
    packages=['lib', 'test'],
    url='',
    license='',
    author='Alexander Graf',
    author_email='graf@genzentrum.lmu.de',
    description='SARS-CoV-2 Sewage quality control',
    scripts=['sewage_quality_control.py'],
    install_requires=[
        'matplotlib == 3.4.3',
        'numpy == 1.23.5',
        'pandas == 2.0.2',
        'pyarrow == 12.0.1',
        'python_dateutil == 2.8.2',
        'scikit_learn == 1.2.0',
        'scipy == 1.9.3',
        'seaborn == 0.12.2',
        'tqdm == 4.65.0'
    ],
    include_package_data=True,
    package_data={'': ['data/*.dat']},
)
