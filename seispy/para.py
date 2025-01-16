from os.path import expanduser, join, exists, dirname
import os
import obspy
from obspy.io.sac import SACTrace
import glob
import configparser
from seispy.io import Query
import re


def _load_station_info(pathname, ref_comp, suffix):
    ex_sac = glob.glob(join(pathname, '*{0}*{1}'.format(ref_comp, suffix)))
    if len(ex_sac) == 0:
        raise FileNotFoundError('no such SAC file in {0}'.format(pathname))
    for sac in ex_sac:
        try:
            ex_tr = SACTrace.read(sac, headonly=True)
            break
        except Exception:
            continue
    if (ex_tr.stla is None or ex_tr.stlo is None):
        raise ValueError('The stlo and stla are not in the SACHeader')
    return ex_tr.knetwk, ex_tr.kstnm, ex_tr.stla, ex_tr.stlo, ex_tr.stel


class StaInfo():
    def __init__(self):
        self.network = ''
        self.station = ''
        self.stla = 0.
        self.stlo = 0.
        self.stel = 0.
        self.channel = '*'
        self.location = '*'
    
    def __repr__(self) -> str:
        # print the station name, channel, location
        return '{0}.{1}.{2}:{}'.format(self.network, self.station, self.channel, self.location)

    def link_server(self, server, user=None, password=None):
        """_summary_

        :param server:  Sever name of FDSN web-service, by default 'IRIS'
        :type server: str
        :param user: username, defaults to None
        :type user: str, optional
        :param password: password, defaults to None
        :type password: str, optional
        """
        try:
            self.query = Query(server, user=user, password=password)
        except Exception as e:
            raise ConnectionError('Cannot connect with {} server'.format(server))

    def get_stainfo(self):
        """Get station information

        :return: Station information
        :rtype: dict
        """
        return self.__dict__

    def load_stainfo(self, pathname, ref_comp, suffix):
        """Load station information from SAC file

        :param pathname: Path to SAC file
        :type pathname: str
        :param ref_comp: Reference component
        :type ref_comp: str
        :param suffix: Suffix of SAC file
        :type suffix: str
        """
        (self.network, self.station, self.stla, self.stlo, self.stel) = _load_station_info(pathname, ref_comp, suffix)

    def get_station_from_ws(self, **kwargs):
        """Get station information from web-service with given network and station or other optional condition.

        
        """
        self.query.get_stations(network=self.network, station=self.station,
                                channel=self.channel, location=self.location, 
                                level='channel', **kwargs)
        if len(self.query.stations) != 1 or len(self.query.stations[0]) != 1:
            raise ValueError('More than one station are fetched. Please set up stricter condition in query')
        self.stla = self.query.stations[0][0].latitude
        self.stlo = self.query.stations[0][0].longitude
        self.stel = self.query.stations[0][0].elevation


class RFPara(object):
    def __init__(self):
        """Parameters for receiver function calculation
        """
        self.datapath = expanduser('~')
        self.rfpath = expanduser('~')
        self.catalogpath = join(dirname(__file__), 'data', 'EventCMT.dat')
        self.pjtpath = 'rfpjt.pkl'
        self.cata_server = 'IRIS'
        self.data_server = 'IRIS'
        self.velmod = 'iasp91'
        self.offset = None
        self.tolerance = 210
        self.dateformat = '%Y.%j.%H.%M.%S'
        self.date_begin = obspy.UTCDateTime('19760101')
        self.date_end = obspy.UTCDateTime.now()
        self.magmin = 5.5
        self.magmax = 10
        self.depthmin = 0.
        self.depthmax = 800.
        self.dismin = 30
        self.dismax = 90
        self.ref_comp = 'BHZ'
        self.suffix = 'SAC'
        self.noisegate = 5
        self.noiselen = 50
        self.phase = 'P'
        self.gauss = 2.0
        self.target_dt = 0.01
        self.time_before = 10
        self.time_after = 120
        self.freqmin = 0.05
        self.freqmax = 1
        self.comp = 'RTZ'
        self.decon_method = 'iter'
        self.wlevel = 0.05
        self.itmax = 400
        self.minderr = 0.001
        self.criterion = None
        self.only_r = False
        self.rmsgate = None
        self.switchEN=False
        self.reverseE=False
        self.reverseN=False
        self.use_remote_data=False
        self.data_server_user = None
        self.data_server_password = None
        self.n_proc = 1
        self.stainfo = StaInfo()

    def get_para(self):
        return self.__dict__

    def __str__(self):
        head = ['{}: {}'.format(k, v) for k, v in self.__dict__.items()]
        return '\n'.join(head)

    def _check_date_range(self):
        if self.date_begin > self.stainfo.query.stations[0][0].end_date or \
            self.date_end < self.stainfo.query.stations[0][0].start_date:
            raise ValueError('The date range set does not overlap with'
                             'the date range of the fetched data')
        if self.date_begin < self.stainfo.query.stations[0][0].start_date:
            self.date_begin = self.stainfo.query.stations[0][0].start_date
        if self.date_end > self.stainfo.query.stations[0][0].end_date:
            self.date_end = self.stainfo.query.stations[0][0].end_date

    @property
    def datapath(self):
        return self._datapath

    @datapath.setter
    def datapath(self, value):
        if not isinstance(value, str):
            raise TypeError('datapath should be \'str\' type not \'{0}\''.format(type(value)))
        else:
            self._datapath = value

    # Temporary function. It will be removed after 2~3 versions
    @property
    def server(self):
        return self.cata_server
    
    @server.setter
    def server(self, value):
        self.cata_server = value

    @property
    def rfpath(self):
        return self._rfpath

    @rfpath.setter
    def rfpath(self, value):
        if not isinstance(value, str):
            raise TypeError('rfpath should be \'str\' type not \'{0}\''.format(type(value)))
        elif not exists(value):
            try:
                os.makedirs(value)
            except Exception as e:
                Warning('Cannot create rfpath of {0}\n with error: {1}'.format(value, e))
            finally:
                self._rfpath = value
        else:
            self._rfpath = value

    @property
    def catalogpath(self):
        return self._catalogpath

    @catalogpath.setter
    def catalogpath(self, value):
        if not isinstance(value, str):
            raise TypeError('catalogpath should be \'str\' type not \'{0}\''.format(type(value)))
        self._catalogpath = value

    @property
    def criterion(self):
        return self._criterion

    @criterion.setter
    def criterion(self, value):
        if self.phase == 'S':
            if value == '' or value is None:
                self._criterion = None
            elif value.lower() == 'lab':
                self._criterion = 'lab'
            else:
                raise ValueError('criterion should be \'lab\'')
        else:
            if ''.join(sorted(self.comp.upper())) == 'LQT':
                self._criterion = None
            elif value == '' or value is None:
                self._criterion = None
            elif value.lower() not in ['crust', 'mtz']:
                raise ValueError('criterion should be string in \'crust\' or  \'mtz\'')
            else:
                self._criterion = value.lower()

    @property
    def decon_method(self):
        return self._decon_method
    
    @decon_method.setter
    def decon_method(self, value):
        if value not in ['iter', 'water']:
            raise ValueError('decon_method must be in \'iter\' or \'water\'')
        else:
            self._decon_method = value

    @classmethod
    def read_para(cls, cfg_file):
        """Read parameters from configure file
        :param cfg_file: Path to configure file
        :type cfg_file: str
        """
        cf = configparser.RawConfigParser(allow_no_value=True)
        pa = cls()
        try:
            cf.read(cfg_file)
        except Exception:
            raise FileNotFoundError('Cannot open configure file %s' % cfg_file)
        sections = cf.sections()
        for key, value in cf.items('path'):
            if value == '':
                continue
            elif key == 'datapath':
                pa.datapath = value
            elif key == 'rfpath':
                pa.rfpath = value
            elif key == 'catalogpath':
                pa.catalogpath = value
            else:
                exec('pa.{} = value'.format(key))
        sections.remove('path')
        if 'fetch' in sections:
            for key, value in cf.items('fetch'):
                if key == 'use_remote_data':
                    pa.use_remote_data = cf.getboolean('fetch', 'use_remote_data')
                elif key == 'data_server_user':
                    if value == '':
                        continue
                    else:
                        pa.data_server_user = value
                elif key == 'data_server_password':
                    if value == '':
                        continue
                    else:
                        pa.data_server_password = value
                elif key == 'data_server':
                    pa.data_server = value
                elif key == 'n_proc':
                    pa.n_proc = cf.getint('fetch', 'n_proc')
                else:
                    exec('pa.stainfo.{} = value'.format(key))
            sections.remove('fetch')
        for sec in sections:
            for key, value in cf.items(sec):
                if key == 'date_begin':
                    pa.date_begin = obspy.UTCDateTime(value)
                elif key == 'date_end':
                    pa.date_end = obspy.UTCDateTime(value)
                elif key == 'offset':
                    try:
                        pa.offset = cf.getfloat(sec, 'offset')
                    except:
                        pa.offset = None
                elif key == 'itmax':
                    pa.itmax = int(value)
                elif key == 'only_r':
                    pa.only_r = cf.getboolean(sec, 'only_r')
                elif key == 'criterion':
                    pa.criterion = value
                elif key == 'decon_method':
                    pa.decon_method = value
                elif key == 'gauss':
                    try:
                        pa.gauss = cf.getfloat(sec, 'gauss')
                    except:
                        pa.gauss = [float(v) for v in re.split(' |,', value)]
                elif key == 'rmsgate':
                    try:
                        pa.rmsgate = cf.getfloat(sec, 'rmsgate')
                    except:
                        pa.rmsgate = None
                else:
                    try:
                        exec('pa.{} = {}'.format(key, float(value)))
                    except ValueError:
                        exec('pa.{} = value'.format(key))
        return pa