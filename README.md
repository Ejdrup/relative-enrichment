# Relative Enrichment
By Aske Lykke Ejdrup
<br>
Last updated on 2022/02/07.

Repository method presented in the paper:
<br>
**Relative enrichment – a density-based colocalization measure for single-molecule localization microscopy**

## Function description
The relative enrichment package consists of two main modules: one to compute the relative enrichment of each reference localization (RE module), and one to bin these values (bin_RE module).
<br><br>
Each of these modules have two main functions, with another two supporting functions for the RE module. All of them are stored in the RE_function.py file, which can be loaded as a package when stored in the current directory by calling:
`import RE_function as re`.

### RE module
The RE module computes the relative enrichment of each reference region based on the two datasets passed to the function. The 2D function is called RE, and can be run by typing `re.RE()`, whereas the 3D version is called RE3D, and can be executed with `re.RE3D()`.
<br><br>
The `re.RE()` takes three inputs:
1. ch1_points - The reference species data set (Nx2).
2. ch2_points - The primary species data set (Mx2).
3. verbose - Whether to output progress (`True` or `False`).

The `re.RE3D()` takes similar inputs, but with Nx3 and Mx3 dimensions.
<br><br>
Both functions output four elements:
<br>
1. n_point_ch2_in_region - Number of primary species in each reference region.
2. sorted_region_area (2D)/sorted_region_volume (3D) - Area or volume of each reference region.
3. first_order_mean_distance - Mean distance to all first order neighbors.
4. bool_index - Boolean area specifying the reference regions with finite size.

The last element is always the same length as the input reference dataset, whereas reference regions with nonfinite size are excluded in the first three outputs. This element is useful for color-coded qualitative plotting.

### bin_RE module
The bin_RE module bins the values outputted by the RE module by either nearest neighbor distance (`bin_RE`) or area/volume of the reference region (`bin_RE_area`).
<br><br>


## Example: simulated vesicle data
Description 