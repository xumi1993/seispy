# -*- coding:utf-8 --

import obspy
import numpy
from scipy.signal import detrend


class rmrespf:
    '''
    remove resp of seismic signal
    '''
    def __init__(self, st, respf, a1, a2, b1, b2):
        '''
        <st> is signal in obspy parttern
        <respf> is the name of respfile
        '''
        self.st    = st
        self.respf = respf
        self.pre_filt = (a1, a2, b1, b2)
        return
#        <a1><a2><b1><b2> are the frequency boundary of the filter

    def rmsimulate(self, units = 'DIS'):
        '''
        <units>:  Units to return response in ('DIS', 'VEL' or ACC)

        '''
        self.st = remt(self.st)
        self.seedresp = {'filename': self.respf,  # RESP filename
                    # when using Trace/Stream.simulate() the "date" parameter can
                    # also be omitted, and the starttime of the trace is then used.

                    #'date': date,

                    # Units to return response in ('DIS', 'VEL' or ACC)
                    'units': units
                   }
        # Remove instrument response using the information from the given RESP file

        self.st.simulate(paz_remove=None, pre_filt=self.pre_filt, seedresp=self.seedresp)
        return self.st

