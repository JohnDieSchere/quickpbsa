# quickpbsa
## Fast and Complete Photobleaching Step Analysis

**Author:** Johan Hummert  
**Organization:** Herten Lab for Single Molecule Spectroscopy, University of Birmingham, UK  
**License:** GPLv3  
**Version:** 2020.0.1

Python Package providing a framework for photo-bleaching step analysis. The details of the algorithm used and extensive validation with experimental data are described in a bioRxiv preprint ().
please cite this publication if you found this package useful.

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
```sh
pip install quickpbsa
```
Alternatively you can clone / download the git repository and place the directory quickpbsa/quickpbsa in your $PYTHONPATH.

## Getting started

If you already have photobleaching traces which you would like to analize, running the analysis is a one-liner:
```python
import quickpbsa as pbsa
pbsa.pbsa_file(file, threshold, maxiter)
```
For this to work the `file` should be a .csv file, where each row is one photobleaching trace, which ideally, but not necessarily, should be background subtracted. Then there are two additional parameters to set in the analysis, `threshold` and `maxiter`. There are many additional optional parameters to set, detailed below, but in many cases the default parameters should be alright.

`threshold` should be set to approximately half the intensity difference of a typical photobleaching step. This is most easily accomplished by plotting a few traces and finding steps towards the end of the trace, where steps can most easily be found by eye.

`maxiter` is the maximum number of iterations (i.e. maximum number of steps found). This should be significantly higher than the expected number of steps. In the validation experiments we performed this would typically be set to 200 for samples with up to 35 fluorophores.

### Structure of the result

The result is exported as a _result.csv file with 7 rows per photobleaching trace. The final result of the complete analysis for each trace is the one with the type 'fluors_full' in the output file. Furthermore traces where the column 'flag' is not 1 should be discarded. The complete structure of the output file and the individual flags are detailed in the Concept below. Examples on how to use the output are also provided under Examples.

### with trace extraction

If you have a stack of .tif images and want to extract and analyze traces from it make sure that you have [tifffile](https://pypi.org/project/tifffile/) installed. Then there are two options to get traces from your image stack, both of which will yield a _difference.csv file containing background corrected traces ready for the analysis:

#### based on localization

If the structures from which you want to extract fluorophore numbers are diffraction-limited, you can extract traces based on a .csv file with coordinates (in nanometers) for each diffraction limited spot. The file should contain at least 2 columns named 'x \[nm\]' and 'y \[nm\]'. Our recommended way of obtaining such a file is the [Fiji](https://fiji.sc/) plugin [ThunderSTORM](https://github.com/zitmen/thunderstorm).

The trace extraction can then be accomplished with
```python
import quickpbsa as pbsa
pbsa.trace_extraction.extract_traces_localization(tiffstack, locfile, r_peak, r_bg1, r_bg2, min_dist)
```
with `tiffstack` being the path to your image stack and `locfile` the path to the localization file. `r_peak` is the radius (in pixels) of the area around the localization from which the trace is extracted. `r_bg1` and `r_bg2` define a ring around the localization from which the background for background correction is extracted. `min_dist` is the minimum distance from one localization to the next. Localizations which are spaced less than `min_dist` apart are not considered in the trace extraction.

#### based on a selection mask

If you have larger structures from which to extract photobleaching traces and fluorophore numbers, you can use a mask image which should be an 8bit Tiff with a white selection on black background. The traces are then extracted from non-connected white regions of interest (ROIs or RsOI ...?) in the mask image.
```python
import quickpbsa as pbsa
pbsa.trace_extraction.extract_traces_mask(tiffstack, maskfile, dist, r_bg)
```
where a background ROI with a distance of `dist` (in pixels) to the selected ROI and a width of `r_bg` is defined for background correction.

## Concept

This photobleaching step analysis is a combination of a preliminary step detection and a following refinement of the preliminary result based on a bayesian approach. A full run of the analysis after trace extraction consists of 3 parts.

### Preliminary step detection

The preliminary step detection is based on the works of Kalafut and Vischer () ...

The optional parameters of the preliminary step detection can be set by providing a dictionary as an optional argument in `quickpbsa.pbsa_file()`:
```python
import quickpbsa as pbsa
pardict = {'norm': 1000,
           'crop': False}
pbsa.pbsa_file(file, threshold, maxiter, preliminary_optional=pardict)
```
Possible parameters are:

|    Parameter      |      Default      | Explanation |
|----------|:-------------:|------:|
| `norm` |  1 | Trace intensity value are divided by norm, mainly for visualization |


### Filtering of traces

The parameters of filtering can be specified in

Based on the result of the preliminary step detection traces are excluded from the analysis. Assuming that the last two steps are correctly identified in a majority of traces, traces are 

## Examples

Detailed examples with a more in-depth explanation of the algorithm are available in two jupyter notebooks explaining how the analysis works on an example trace (Examples/Example_Trace.ipynb) and an example tiff-stack (Examples/Example_Stack.ipynb) with the included experimental example data. You will need to install [jupyter](https://jupyter.org/) and [matplotlib](https://matplotlib.org/) to run the examples.

## Function Reference





