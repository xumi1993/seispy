"""
Created to replace class RFStation
containing 2 class:

1. Trace_Clust: similar to obspy.core.Stream. but traces are in SACTrace format, and list is at self.data
2. RFSta: rfs of same station. Created to replace class RFStation. use self.check_unit to do header checks


Notes:
    SACTrace can be convert to obspy.core.Trace object or convert from Trace
    compare to RFStation, RFSta add list-like operation, and discard phase option
    phase or datar/datat options are discarded
"""

from collections import UserList
from os.path import isfile, exists

import numpy as np

import obspy
from obspy.io.sac import SACTrace
from obspy.core import Stream


class Trace_Clust(UserList):
    """
    user can accese traces from self.data or like a default list

    .. rubric:: Example
    
    >>> st= obspy.read();tt=Trace_Clust()
    >>> for _i in st:
    ...   _i.stats.sac ={};
    ...   _i.stats.sac["stla"]=10.0
    ...   _i.stats.sac["stlo"]=20.0
    >>> st = [SACTrace.from_obspy_trace(_i) for _i in st]
    >>> len(st)
    3

    you can directly append trace or traces in SACTrace format
    >>> tt.append(st[0]);len(tt)
    1
    >>> tt.append(st);len(tt)
    4

    [] is available either
    >>> tt[0].npts
    3000

    for iter is also available
    >>> [_i.npts for _i in tt]
    [3000, 3000, 3000, 3000]
    """

    def __init__(self, input_list=[]):
        """
        init self.data, from list
        """
        if isinstance(input_list, Trace_Clust):
            self.data=input_list.data
        elif isinstance(input_list, list) or isinstance(input_list, SACTrace):
            self.data=[]
            self.append(input_list)
        else:
            self.data=[]

    def append(self, item):
        if isinstance(item, SACTrace):
            self.data.append(item)
        elif isinstance(item, list) or isinstance(item, UserList):
            for _i in item:
                self.append(_i)

    def __add__(self, other):
        if isinstance(other, Trace_Clust):
            return self.__class__(self.data + other.data)
        elif isinstance(other, SACTrace):
            return self.__class__(self.data + list(other))
        return self.__class__(self.data)


class RFSta(Trace_Clust):
    """
    used to replace RFStation for ccp usage,
    contains rfs of same stations, and behave like a default list in python
    many test are show in Trace_Clust sec

    here we only do a set of simple test to show how to read file and
    """
    # this is for compatible between seispy and sac
    _attr_same={
        "stla": "stla", "stlo": "stlo", "stel": "stel",
        "rflength": "npts", "shift": "b", "sampling": "delta",
    }
    # in check_
    # we assume rayp is cal and stored in user0 of sac headers
    _attr_mat={
        "rayp": "user1", "bazi": "baz", "mag": "mag",
        "evla": "evla", "evlo": "evlo", "evdp": "evdp", "dis": "gcarc"
    }
    DEFAULT_PRECISION=0.05

    def from_stream(self, st: Stream):
        """
        import sac traces from obspy.core.Stream object\
        so users can use * or other easier way to import sac traces
        >>> st= obspy.read();tt=Trace_Clust()
        >>> for _i in st:
        ...   _i.stats.sac ={};
        ...   _i.stats.sac["stla"]=10.0
        ...   _i.stats.sac["stlo"]=20.0
        >>> tt=RFSta(); tt.from_stream(st);len(tt)
        3
        """
        self.data=[SACTrace.from_obspy_trace(_i) for _i in st]
        self.check_unit(accor=True)

    def check_unit(self, **kwargs):
        """
        make headers for RFStation
        kwargs marks different branches of check
        contains:
        1. accor ---- check stla,stlo,stel,shift,npts and so on. RF with different shift, delta will be discard
        2. postition ---- check stla,stlo for all station, may merge to 1
        3. gen_mat ---- gen mats contains baz,rayp and so on
        4. accor_seispy ---- gen ev_num and data_prime

        replace read_sample
        >>> st= obspy.read();tt=Trace_Clust()
        >>> for _i in st:
        ...   _i.stats.sac ={};
        ...   _i.stats.sac["stla"]=10.0
        ...   _i.stats.sac["stlo"]=20.0
        >>> tt=RFSta(); tt.from_stream(st);tt.check_unit(accor=True,postion=True, gen_mat=True, accor_seispy=True)
        >>> len(tt)
        3
        >>> tt[2].stla = 15.0; tt.check_unit(position = True)
        >>> len(tt)
        2
        """
        if len(self) == 0:
            return -1

        # init attr if not exist
        # _atter_same keeps values that same for all stations
        # like stla stlo, here we use first one in seq to set value
        # one thing for confirm is shift eq -b while cal for rf in seispy
        if kwargs.pop("accor", True):
            self.ev_num=len(self)
            for _i, _j in self._attr_same.items():
                # in case stel is not set in head info
                if not hasattr(self, _i):
                    setattr(self, _i,
                            getattr(self[0], _j, 0.))
            self.shift=-self.shift
            self.staname='{}.{}'.format(
                self[0].knetwk,
                self[0].kstnm
            )
            # npts and shift check
            picks=[]
            for _j, _i in enumerate(self.data):
                if (_i.delta - self[0].delta) > self.DEFAULT_PRECISION \
                        or (_i.b - self[0].b) > 5 * _i.delta:
                    picks.append(_j)
            self.clean(picks)
            picks.clear()

        # use position,shift to roughly constrain
        if kwargs.pop("position", False):
            if not hasattr(self, "stla") or not hasattr(self, "stlo"): self.check_unit("accor")
            precision=kwargs.pop('precision', self.DEFAULT_PRECISION)
            picks=[]
            for _i in range(1, len(self)):
                if abs(self[_i].stla - self.stla) > precision \
                        or abs(self[_i].stlo - self.stlo) > precision \
                        or abs(self[_i].b + self.shift) > 0.2:
                    picks.append(_i)
            self.clean(picks)
            picks.clear()

        # set mat of rayp and bazi
        # not sure list is suitable for this design
        # if self.trace is somehow sorted, there's need to re-arrange this format
        # few value comfirm, dangrous
        if kwargs.pop("gen_mat", False):
            size=len(self)
            for _i, _j in self._attr_mat.items():
                _mat=[]
                for _t in range(0, size):
                    _mat.append(getattr(self[_t], _j))
                setattr(self, _i, _mat)

        # gen matrix for seispy's posible usage
        # warning: This step hasn't set any protection for different sampling rate
        if kwargs.pop("accor_seispy", False):
            length=max([_i.npts for _i in self])
            self.data_prime=np.zeros(
                (len(self), length)
            )
            for _i, _j in enumerate(self.data):
                self.data_prime[_i:_j.npts]=_j.data

    def clean(self, rm_list=[]):
        if not isinstance(rm_list, list): return -1
        tmp=[]
        for _i, _j in enumerate(self.data):
            if _i not in rm_list:
                tmp.append(_j)

        self.data=tmp


if __name__ == "__main__":
    import doctest

    doctest.testmod()
