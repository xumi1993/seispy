# CCP stacking

```python

[FileIO]
    # Path to stations with RF sac files
rfpath = /path/to/RFresult

# Path to station list
stalist = /path/to/sta_all.lst

# Path to the lib of Ps ray parameters. 
# If it's empty the ray parameters of Ps would be assumed as that of P arrival
rayp_lib =

# Output data structure after time to depth
depthdat =  /path/to/RFdepth

# Output data structure after CCP stacking
stackfile = /path/to/stack_data

# Station list used to stack
stack_sta_list = /path/to/stack_sta.lst

# Path to 1D velocity model with 3 columns: depth vp vs
# If it's empty, the IASP91 model will be used in time-to-depth conversion 
velmod =

# Optional, Path to file for searching depth of d410 and d660
peakfile = 

[bin]
# For linear array, wether create bins with a self-adaptive method
adaptive = false

# The shape of bins, circle or rect available 
shape = rect

# Radius of bins in km
# Set to empty for determination of radius with fresnel zone
bin_radius =

# period of S wave (for assuming the radius of fresnel zone)
domperiod = 5

# Width of the profile in km, only works for rectangle bin
width = 100

# sliding or spacing interval of bins in km 
slide_val = 5

[line]
# Coordinate of two end points for ccp_profile
profile_lat1 = 27.5
profile_lon1 = 94
profile_lat2 = 36.5
profile_lon2 = 92

[spacedbins]
# Spaced grid for ccp3d
center_lat = 32
center_lon = 94
half_len_lat = 4
half_len_lon = 4

[depth]
# Max depth for time-to-depth conversion
dep_end = 800
dep_val = 1

[stack]
# Stack RFs from <stack_start> km to <stack_end> km with interval of <stack_val> km
stack_start = 300
stack_end = 750
stack_val = 2

# Samples in bootstrap method
boot_samples = 2000
```