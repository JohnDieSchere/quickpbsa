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

**norm (default 1)** The traces will be divided by `norm` prior to analysis, mainly for visualization. `threshold` should be set lower accordingly.  
**crop (default True)** If `True` the last frames of the trace are cropped for analysis purposes based on `threshold`. Traces are cropped at the last frame where the difference in intensity exceeds half the value of `threshold` + `bgframes`.  
**bgframes (default 500)** How many frames to include in the analysis after the crop point.  
**max_memory (default 2.0)** Maximum available memory for the preliminary step detection. The analysis defaults to a slower, but less memory consuming implementation if the necessary memory exceeds this value. For the default of 4 GB the fast implementation is used for traces with up to ~40000 frames.

### Filtering of traces

Based on the result of the preliminary step detection traces are excluded from the analysis. Assuming that the last two steps are correctly identified in a majority of traces, traces are flagged out. The used flags in the output file are:

| flag | Meaning |
| :--  | :------
| -1   | No steps found in preliminary step detection
| -2   | Background intensity out of bounds
| -3   | Single fluorophore intensity out of bounds
| -4   | Trace goes into negative values
| -5   | Fluorophore number becomes negative
| -6   | interval between final two steps is too short

The optional parameters of the preliminary step detection can be set by providing a dictionary as an optional argument in `quickpbsa.pbsa_file()`:
```python
import quickpbsa as pbsa
pardict = {'subtracted': False,
           'percentile_step': [20, 80]}
pbsa.pbsa_file(file, threshold, maxiter, filter_optional=pardict)
```
Possible parameters are:

**subtracted (default True)** If `True` it is assumed that traces are background corrected. This sets the bounds on the background intensity to `[-threshold, threshold]` and the default lower bound on the single fluorophore intensity to `threshold`. If `False` the bounds on the background intensity are set based on the minimum background intensity in the dataset `min_bg`: `[min_bg, min_bg + threshold]`. If `False` the default lower bound on the single fluorophore intensity is also `min_bg + threshold`.  
**percentile_step (default 90)** Sets the bounds on the single fluorophore intensity. If one value is provided, the upper bound on the single fluorophore intensity is set at this percentile. If two values are provided, as in the example above, lower and upper bounds are set at the percentiles respectively.  
**length_laststep (default 20)** Minimum number of frames between the last two steps.

### Step refinement

The step refinement is based on the posterior as defined in ...

Most of the optional parameters aim to reduce runtime by reducing the number of possible step arrangements to test. The optional parameters of the step refinement can also be set by providing a dictionary as an optional argument in `quickpbsa.pbsa_file()`:
```python
import quickpbsa as pbsa
pardict = {'mult_threshold': 1.5,
           'combcutoff': int(5e6)}
pbsa.pbsa_file(file, threshold, maxiter, refinement_optional=pardict)
```
Possible parameters are:

**multstep_fraction (default 0.5)** Maximum fraction of steps with an occupancy higher than 1.  
**nonegatives (default False)** If `True`, no negative double steps are considered. This means that arrangements where 2 or more fluorophore turn back on at the same time are not considered.  
**mult_threshold (1.0)** Only steps where the difference in the mean is above `mult_threshold` multiplied by the last fluorophore intensity are considered as steps with occupancy higher than 1.  
**combcutoff (default 2000000)** Maximum number of arrangements to test. If this is exceeded, the trace is flagged out with flag `-7`. If this happens a lot, consider increasing this value, which will increase runtime.  
**splitcomb (default 30000)** How many arrangements to test simultaneously (vectorized). On systems with a large memory this can be increased to speed up the analysis.  
**maxmult (default 5)** Maximum considered occupancy, i.e. how many fluorophores can bleach simultaneously.  
**maxadded (default 10)** Maximum number of added single steps if no steps are removed to yield an improved posterior.  
**lambda (default 0.1)** Hyperparameter $\lambda$ in equation (1).  
**gamma0 (default 0.5)** Hyperparameter $\gamma_0$ in equation (1).

## Examples

Detailed examples with a more in-depth explanation of the algorithm are available in two jupyter notebooks explaining how the analysis works on an example trace (Examples/Example_Trace.ipynb) and an example tiff-stack (Examples/Example_Stack.ipynb) with the included experimental example data. You will need to install [jupyter](https://jupyter.org/) and [matplotlib](https://matplotlib.org/) to run the examples.





