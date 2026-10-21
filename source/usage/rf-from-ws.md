# Fetch seismic data from web-service and calculate RF automatically

After Seispy v1.3.0, users can calculate RFs with a specified network and station by fetching data from an [FDSN web service](https://www.fdsn.org/webservices/). This section shows an example of calculating PRFs after fetching station and event information from a web service.

```{note}
This notebook can be downloaded as **{download}`rf-from-ws.ipynb <../_static/rf-from-ws.ipynb>`**
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

Before running this script, we can visually search stations from web-service portals, such as [GFZ webdc3](http://eida.gfz-potsdam.de/webdc3/) or [IRIS GMap](https://ds.iris.edu/gmap/). The URL of an FDSN web service or shortcut names can be found in [obspy.client.fdsn](https://docs.obspy.org/packages/obspy.clients.fdsn.html). The network name, station name, positions, date range, etc. can be found at these services. Now let's fetch station information using these conditions.

The following example illustrates how to request station information from the Global Seismograph Network (`IU`).

```python
query = Query(server='IRIS') ## Server is the URL of FDSN web service or a shortcut name in the obspy.client.fdsn.
query.get_stations(network='IU', station='U*', level='channel')
#print(query.stations)
```

## Fetch data and calculate RF with different gauss factors

The following example illustrates how to request the `'BH?'` channels, `'00'` location of station Ulaanbaatar (`'ULN'`) of the Global Seismograph Network (`'IU'`) for events between "2013-08-01" and "2013-10-31" (UTC), calculate the RF with 4 gauss factors simultaneously, and save the raw seismic data and RFs.

<h3 id="rfpara">Set the parameters for matching catalog and estimating RF</h3>

All parameters for matching catalog and estimating RF are in the [`RF.para`](#rfpara). These parameters can be set according to user needs. Online catalog (`'cata_server'`) is fetched from the FDSN web service client for ObsPy ([obspy.client.fdsn](https://docs.obspy.org/packages/obspy.clients.fdsn.html)).

```python
rf = RF()
rf.para.data_server = 'IRIS'
rf.para.cata_server = 'IRIS'
rf.para.stainfo.network = 'IU'
rf.para.stainfo.station = 'ULN'
rf.para.stainfo.channel = 'BH?'
rf.para.stainfo.location = '00'
rf.para.datapath = './Data/{}.{}'.format(rf.para.stainfo.network, rf.para.stainfo.station)
rf.para.use_remote_data = True
rf.para.ref_comp ='BHZ'
rf.para.phase = 'P'
rf.para.noisegate = 1
rf.para.magmin = 5.8
rf.para.gauss = [0.5, 1.0, 1.5, 2.0] ##RF with different Gauss factor will be calculated simultaneously.
rf.para.rmsgate = 0.4
rf.para.freqmin = 0.05
rf.para.freqmax = 2.0
rf.para.comp = 'RTZ'
rf.para.date_begin = UTCDateTime('20130801')
rf.para.date_end = UTCDateTime('20131031')
```

### Load station information and search events

- Fetch the station information from the data server ([`data_server`](#rfpara)).
- Search the event information from the catalog server ([`cata_server`](#rfpara)).

Here, we use the `'IRIS'` client. Available catalogs in IRIS are listed in [IRIS DMC FDSNWS event Web Service](https://service.iris.edu/fdsnws/event/1/catalogs), such as `'ISC'`, `'NEIC PDE'` and `'GCMT'`.

```python
rf.load_stainfo()
rf.search_eq(catalog='NEIC PDE')
#print(rf.eq_lst) ##The matched event lists are listed.
```

### Match catalog and fetch seismic data

Match events and fetch seismic data with parameters such as the data type (`'SAC'`) and date format `'%Y.%j.%H.%M.%S'` set in the [`RF.para`](#rfpara).

```python
rf.match_eq()
```

### Calculate RF

- Remove the linear trend (`detrend`) and apply a bandpass filter (`filter`) to the data. The frequencies for the bandpass filter are set in the [`RF.para`](#rfpara) ([`'para.freqmin'`](#rfpara) and [`'para.freqmin'`](#rfpara)).
- Mark phase arrivals with the TauP server; the velocity model (`'para.velmod'`) can be set in the [`RF.para`](#rfpara).
- Rotate the seismic data to `'RTZ'` or `'LQT'` and delete events with an SNR lower than [`'para.noisegate'`](#rfpara).
- Save the raw SAC data downloaded from the web server.
- Trim the RF between [`para.time_before`](#rfpara) and [`para.time_after`](#rfpara).
- Deconvolve to obtain RFs with different gauss factors. Available deconvolution methods (`para.decon_method`) are time-domain iterative deconvolution (`'iter'`) and frequency-domain water-level deconvolution (`'water'`).

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

Save the RFs calculated with different Gauss factors.

```python
for ff in rf.para.gauss:
    rf.para.rfpath = './RFresult/F{:.1f}/{}.{}'.format(
        ff, rf.para.stainfo.network, rf.para.stainfo.station)
    rf.saverf(ff)
```
