from setuptools import setup

setup(
   name='blenderspike_py',
   version='1.0',
   description='A companion module for BlenderSpike to export NEURON simulations into pickle',
   author='Artem Kirsanov',
   author_email='ArtemKirsanov2606@gmail.com',
   packages=['blenderspike_py'], 
   install_requires=['numpy', 'scipy'],
)