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
If you would stack PRFs in a great depth (e.g., D410 or D660), the raypara of P arrival cannot represent that of Ps phase. Therefore, 
a library file of Ps ray parameters in different depths and epicentral distance is required in calculating Ps-P time difference.
