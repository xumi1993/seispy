from seispy.rf import RF, datestr2regex, SACFileNotFoundError
from seispy.eq import EQ
import numpy as np
from obspy import UTCDateTime
from datetime import timedelta
import pandas as pd
import glob
from os.path import join
import re
import sys
import obspy


def match_eq(eq_lst, pathname, logger, ref_comp='Z', suffix='SAC', offset=0,
             tolerance=1, dateformat='%Y.%j.%H.%M.%S'):
    pattern = datestr2regex(dateformat)
    ref_eqs = glob.glob(join(pathname, '*{0}*{1}'.format(ref_comp, suffix)))
    if len(ref_eqs) == 0:
        raise SACFileNotFoundError(join(pathname, '*{0}*{1}'.format(ref_comp, suffix)))
    sac_files = []
    for ref_sac in ref_eqs:
        try:
            datestr = re.findall(pattern, ref_sac)[0]
        except IndexError:
            raise IndexError('Error data format of {} in {}'.format(dateformat, ref_sac))
        if isinstance(offset, (int, float)):
            sac_files.append([datestr, UTCDateTime.strptime(datestr, dateformat), -offset])
        elif offset is None:
            try:
                tr = obspy.read(ref_sac)[0]
            except TypeError:
                continue
            sac_files.append([datestr, tr.stats.starttime-tr.stats.sac.b, tr.stats.sac.o])
        else:
            raise TypeError('offset should be int or float type')
    new_col = ['data', 'datestr']
    eq_match = pd.DataFrame(columns=new_col)
    for datestr, b_time, offs in sac_files:
        date_range_begin = b_time + timedelta(seconds=offs - tolerance)
        date_range_end = b_time + timedelta(seconds=offs + tolerance)
        results = eq_lst[(eq_lst.date > date_range_begin) & (eq_lst.date < date_range_end)]
        if len(results) != 1:
            continue
        try:
            this_eq = EQ(pathname, datestr, suffix)
        except Exception as e:
            logger.RFlog.error(''.format(e))
            continue
        this_eq.get_time_offset(results.iloc[0]['date'])
        this_df = pd.DataFrame([[this_eq, datestr]], columns=new_col, index=results.index.values)
        eq_match = pd.concat([eq_match, this_df])
    ind = eq_match.index.drop_duplicates(keep=False)
    eq_match = eq_match.loc[ind]
    return pd.concat([eq_lst, eq_match], axis=1, join='inner')


class ReRF(RF):
    def __init__(self, finallist, cfg_file=None, log=None):
        super().__init__(cfg_file, log)
        self.eq_lst = self.read_finallist(finallist)
    
    def read_finallist(self, finallist):
        """Read final list from RF path.

        Parameters
        ----------
        finallist : str
            Path to final list.

        Returns
        -------
        pandas.DataFrame
            Data frame of earthquake information.
        """
        self.logger.RFlog.info('Read event info from {}'.format(finallist))
        dtype = {'names': ('evt', 'phase', 'evlat', 'evlon', 'evdep', 'dis', 'bazi', 'rayp', 'mag', 'f0'),
                 'formats': ('U20', 'U20', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')}
        self.event, self.phase, self.evla, self.evlo, self.evdp, self.dis, self.bazi, self.rayp, self.mag, self.f0 = \
            np.loadtxt(finallist, dtype=dtype, unpack=True, ndmin=1)
        cols = ['date', 'evla', 'evlo', 'evdp', 'dis', 'bazi', 'rayp', 'mag', 'magtype']
        lst = []
        for i, evt_str in enumerate(self.event):
            evtime = UTCDateTime.strptime(evt_str, '%Y.%j.%H.%M.%S')
            lst.append([evtime, self.evla[i], self.evlo[i], self.evdp[i],
                       self.dis[i], self.bazi[i], self.rayp[i], self.mag[i], 'mw'])
        return pd.DataFrame(lst, columns=cols)
    
    def match_eq(self, **kwarg):
        try:
            self.logger.RFlog.info('Match SAC files')
            self.eqs = match_eq(self.eq_lst, self.para.datapath, self.logger,
                                ref_comp=self.para.ref_comp, suffix=self.para.suffix,
                                offset=self.para.offset, tolerance=self.para.tolerance,
                                dateformat=self.para.dateformat)
        except Exception as e:
            self.logger.RFlog.error('{0}'.format(e))
            raise e
        if self.eqs.shape[0] == 0:
            self.logger.RFlog.warning('No earthquakes matched, please check configurations.'.format(self.eqs.shape[0]))
            sys.exit(1)
        else:
            self.logger.RFlog.info('{0} earthquakes are matched'.format(self.eqs.shape[0]))

    def write_list(self):
        path = join(self.para.rfpath, '{}.{}finallist.dat'.format(self.stainfo.network, self.stainfo.station))
        self.logger.RFlog.info('Writting event info to {}'.format(path))
        with open(path, 'w') as f:
            for i, row in self.eqs.iterrows():
                f.write('{} {} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.5f} {:.3f} {:.3f}\n'.format(
                    row['date'].strftime('%Y.%j.%H.%M.%S'), self.para.phase, row['evla'], row['evlo'],
                    row['evdp'], row['dis'], row['bazi'], row['rayp'], row['mag'], self.para.gauss
                ))
