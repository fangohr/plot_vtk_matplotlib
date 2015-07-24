from setuptools import setup

setup(
    version='0.1',
    name="plot_vtk_matplotlib",
    packages=['plot_vtk_matplotlib'],
    license='BSD',
    description='Library to plot VTK files in 2D with Matplotlib',
    # We comment this to avoid upgrading matplotlib when
    # instaling the library. Could be used in the future.
    # install_requires=['scipy', 'matplotlib'],
    author_email='d.i.cortes@soton.ac.uk'
)
