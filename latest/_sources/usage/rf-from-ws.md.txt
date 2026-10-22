# Fetch seismic data from web-service and calculate RF automatically

After Seispy v1.3.0, user can calculate RFs with specified network and station which can be fetched from [FDSN web service](https://www.fdsn.org/webservices/). This section shows a example to calculate PRFs with fetching station and event information from web service.

```{note}
This notebook can be downloaded as **{nb-download}`rf-from-ws.ipynb`**
```


```python
import os
from seispy.rf import RF
from seispy.io import Query
from obspy import UTCDateTime
#import pytest
import glob
```

## Get information of stations
Before running this script, we can visually search stations from portal of web-service, such as [GFZ webdc3](http://eida.gfz-potsdam.de/webdc3/) or [IRIS GMap](https://ds.iris.edu/gmap/). The URL of FDSN web service or shortcut names can be found in [obspy.client.fdsn](https://docs.obspy.org/packages/obspy.clients.fdsn.html). The network name, station name, positions, date range, etc. can be found at these services. Now let's fetch station information using these conditions.  
   
The following example illustrates how to request the station information from the Global Seismograph Network("`IU`").


```python
query = Query(server='EARTHSCOPE') ## Server is the URL of FDSN web service or a shortcut name in the obspy.client.fdsn.
query.get_stations(network='IU', station='U*', level='channel')
#print(query.stations)
```

## Fetch data and calculate RF with different gauss factors
    
   The following example illustrates how to request the `'BH?'` channels, `'00'` location of station Ulaanbaatar (`'ULN'`) of the Global Seismograph Network(`'IU'`) for events between "2013-08-01"  and "2013-10-31" (UTC), calculate the RF with 4 gauss factors simultaneously, and save the raw seismic data and RFs.   
      

<h3 id="rfpara">Set the parameters for matching catalog and estimating RF</h3>  
   
   All parameters for matching catalog and estimating RF are in the [`RF.para`](#rfpara).  These parameters can be set according to user needs. 
   Online catalog (`'cata_server'`) is fetched from the FDSN web service client for ObsPy ([obspy.client.fdsn](https://docs.obspy.org/packages/obspy.clients.fdsn.html)). 


```python
rf = RF()
rf.para.data_server = 'EARTHSCOPE'
rf.para.cata_server = 'USGS'
rf.para.stainfo.network = 'IU'
rf.para.stainfo.station = 'ULN'
rf.para.stainfo.channel = 'BH?'
rf.para.stainfo.location = '00'
rf.para.datapath = './Data/{}.{}'.format(rf.para.stainfo.network, rf.para.stainfo.station)
rf.para.use_remote_data = True
rf.para.ref_comp ='BHZ'
rf.para.phase = 'P'
rf.para.noisegate = 1
rf.para.n_proc = 4
rf.para.magmin = 5.8
rf.para.gauss = [0.5, 1.0, 1.5, 2.0] ##RF with different Gauss factor will be calculated simultaneously.
rf.para.rmsgate = 0.4
rf.para.freqmin = 0.05
rf.para.freqmax = 2.0
rf.para.comp = 'RTZ'
rf.para.date_begin = UTCDateTime('20130801')
rf.para.date_end = UTCDateTime('20131031')

```

### load station information and search events  
  - Fetch the station information from the data server ([`data_server`](#rfpara))  
  - Search the event information from the catalog server ([`cata_server`](#rfpara)).  
  Here, we use the `'IRIS'` client. Available catalogs in the IRIS are listed in [IRIS DMC FDSNMS event Web Server](https://service.iris.edu/fdsnws/event/1/catalogs), such as `'ISC'`, `'NEIC PDE'` and `'GCMT'`.


```python
rf.load_stainfo()
rf.search_eq()
#print(rf.eq_lst) ##The matched event lists are listed.
```

### Match catalog and fetch seismic data 
  Match events and fetch seismic data with the parameters such as the data type (`'SAC'`)  and dateformat `'%Y.%j.%H.%M.%S'` set in the [`RF.para`](#rfpara).
  


```python
rf.match_eq() 
```

### Calculate RF 
  - Remove the linear trend (`detrend`) and apply a bandpass filter (`filter`) to the data. The frequencies for the bandpass filter are set in the [`RF.para`](#rfpara) ([`'para.freqmin'`](#rfpara) and [`'para.freqmin'`](#rfpara));   
  - Mark phase arrivals with the server of Taup and the velocity mode ([`'para.velmod'`](#rfpara)) can be set in the [`RF.para`](#rfpara);   
  - Rotate the seismic data to `'RTZ'` or`'LQT'` and delete the events with the SNR lower than the [`'para.noisegate'`](#rfpara);   
  - Save the raw SAC data download from the web server;   
  Trim the RF bewteen the times of [`para.time_before`](#rfpara) and [`para.time_after`](#rfpara);    
  - Do deconvolution to obtain the RFs with different gauss factors. Deconvolution methods (`para.decon_method`) of Time-domain iterative deconvolution (`'iter'`) and frequency-domian water-level deconvolution (`'water'`) are available.  


```python
rf.detrend()
rf.filter()
rf.cal_phase()
rf.rotate()
rf.drop_eq_snr()
rf.save_raw_data()
rf.trim()
rf.deconv()
```

### Save the RFs   
  Save the RFs calculating with different Gauss factors. 


```python
for ff in rf.para.gauss:
    rf.para.rfpath = './RFresult/F{:.1f}/{}.{}'.format(
        ff, rf.para.stainfo.network, rf.para.stainfo.station)
    rf.saverf(ff)
```
