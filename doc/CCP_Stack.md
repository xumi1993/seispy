# Stack PRFs by Common Conversion points method

This note will introduce how to stack PRFs by CCP method in two steps.
1. Convert PRFs from time axis to depth axis according to 1D velocity model.
2. search and stack amplitudes in each bin with calculated pierce points falling in

## Preparation
Befor the CCP stacking, we should prepare:
### 1. PRFs data with sac format

- the SACHeader of `b`, `delta` and `npts` are required. please make sure that they are correct. The `b` denote the shift time before P arrival.
- All of SAC files belong to a same station must be put into a same directory with name of the station name.
- A list file with 10 columns is required in this directory, whose format as:
    - `evt`: the datetime part of sacfile name of the PRF. 
    - `phase`: The phase part of sacfile name of the PRF.
    - `evlat`: The latitude of the event.
    - `evlon`: The longitude of the event.
    - `evdep`: The focal depth of the event.
    - `dis`: The epicentral distance between event and station.
    - `bazi`: The back-azimuth from event to station.
    - `rayp`: The ray parameter of the event.
    - `mag`: The magnitude of the event.
    - `f0`: The gauss factor used in computing the PRF.

### 2. (Optional) a lib of Ps ray parameters.
If you would stack PRFs in a great depth (e.g., D410 or D660), the ray parameters of P arrival cannot represent that of Ps phase. Therefore, 
a library file of Ps ray parameters in different depths and epicentral distance is required in calculating Ps-P time difference.
This lib file has specific binary format. You can generate it by command of `gen_rayp_lib`:
```
usage: gen_rayp_lib [-h] -d DIS_STR -e DEP_STR [-l LAY_STR] [-o OUT_PATH]

Gen a ray parameter lib for Pds phases

optional arguments:
  -h, --help   show this help message and exit
  -d DIS_STR   Distance range as 'min_dis/max_dis/interval'
  -e DEP_STR   Event depth range as 'min_dep/max_dep/interval'
  -l LAY_STR   layers range as 'min_layer/man_layer'
  -o OUT_PATH  Out path to Pds ray parameter lib
```