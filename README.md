# quickpbsa

Python Package providing a framework for photo-bleaching step analysis. The details of the algorithm used and extensive validation with experimental data are described in a bioRxiv preprint ().
please cite this publication if you found this package useful.

**Author:** Johan Hummert  
**Organization:** Herten Lab for Single Molecule Spectroscopy, University of Birmingham, UK  
**License:** GPLv3  
**Version:** 2020.0.1

## Dependencies

Although the package was tested with specific versions of these packages, other relatively new versions will likely work as well. If you have issues with newer versions of these packages get in touch. If you have issues with older versions please consider updating.

- **Python 3:**
	Tested with python 3.8
- **numpy:**
	Tested with numpy-1.18.2
- **scipy:**
	Tested with scipy-1.4.1
- **pandas:**
	Tested with numpy-1.0.3
- **sympy:**
	Tested with 1.5.1
	
### Optional dependendencies for trace extraction

If you want to use the package not only for analysis of photobleaching traces but also for trace extraction from .tiff stacks there are additional dependencies:

- **tifffile:**
	Tested with tifffile-2019.2.10
- **matplotlib:**
	Tested with matplotlib-3.2.1

## Installation

The recommended way to install is via pip:
```python
pip install quickpbsa
```
Alternatively you can clone / download the git repository and place the directory quickpbsa/quickpbsa in your $PYTHONPATH.

## Getting started

If you already have photobleaching traces which you would like to analize, running the analysis is a one-liner:
```python
import quickpbsa as pbsa
pbsa.pbsa_file(file, threshold, maxiter)
```
For this to work the ```file``` should be a .csv file, where each row is one photobleaching trace, which ideally, but not necessarily, should be background subtracted. Then there are two additional parameters to set in the analysis, ```threshold``` and ```maxiter```.

```threshold``` should be set to approximately half the intensity difference of a typical photobleaching step. This is most easily accomplished by plotting a few traces and finding steps towards the end of the trace, where steps can most easily be found by eye.

```maxiter``` is the maximum number of iterations (i.e. maximum number of steps found). This should be significantly higher than the expected number of steps. In the validation experiments we performed this would typically be set to 200 for samples with up to 35 fluorophores.

### 

## Concept


## Examples

Detailed examples with a more in-depth explanation of the algorithm are available in two jupyter notebooks explaining how the analysis works on an example Trace (Examples/Example_Trace.ipynb) and an example tiff-stack (Examples/Example_Stack.ipynb) with the included (experimental) example data. You will need to install jupyter-notebook (https://jupyter.org/) and matplotlib (https://matplotlib.org/) to run those examples.

## Function Reference




