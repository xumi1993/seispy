# RF calculation

```
[path]
# Data path with event records
datapath = Data.CB.NJ2

# Output RF path
rfpath = ./RFresult/CB.NJ2

# Catalog path
catalogpath = events_all.dat

path to project file
pjtpath = rfpjt.pkl

[search_eq]
# Date range
date_begin = 20140101
date_end = 20150723

# Catalog server, defaults to IRIS
catalog_server = IRIS

# Events limitation
magmin = 5.5
magmax = 10
dismin = 30
dismax = 90

[match_eq]
# Date format in SAC filenames
dateformat = %Y.%j.%H%M%S

# Reference component in SAC filenames
ref_comp = .BHZ.

# Suffix of SAC filenames
suffix = SAC

# Relative time and tolerance
offset = 0
tolerance = 10

[snr]
# Threshold of SNR
noisegate = 1

# Time length of SNR calculation
noiselen = 50

[filter]
# Low and high cut off frequency
freqmin = 0.05
freqmax = 2

[rotation]
# Target components, TRZ or LQT available
comp = TRZ

[trim]
# Time before and after P arrival for trim
time_before = 10
time_after = 120

[decon]
# Deconvolution method, iter or water available
decon_method = iter

# Gauss factor
gauss = 2

# Iteration times and min error
itmax = 400
minderr = 0.001

# Water level
wlevel = 0.05

[save]
# DT for resample
target_dt = 0.01

# Calculate only R component, yes or no available
only_r = no

# Criterion for saving RFs, non-set, crust or mtz available
criterion = crust

# Max RMS for RFs
rmsgate = 0.25
```