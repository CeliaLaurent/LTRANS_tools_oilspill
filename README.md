# LTRANS_tools python library to post-proc oilspill simulations

- to use this library define first the ambient variable:
`export LTRANS_tools=/path/to/this/folder/` 
- specify paths to input files in `pyscripts/compute_oilspill.py` 
- run the script `pyscripts/compute_oilspill.py` with python3.
-This script produces a binary file whose name and path must be specified in `pyscripts/plot_oilspill_hazardmap.py`
- run the script `pyscripts/plot_oilspill_hazardmap.py` with python3 to produce the hazard maps.
