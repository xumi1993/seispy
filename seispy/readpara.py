from configparser import ConfigParser, RawConfigParser
from seispy.rf import para, rf
import obspy


def readpara(cfg_file):
    pa = para()
    cf = RawConfigParser()
    cf.read(cfg_file)
    for key, value in cf.items('path') + cf.items('para'):
        if key == 'date_begin' or key == 'date_end':
            pa.__dict__[key] = obspy.UTCDateTime(value)
        elif key == 'only_r':
            pa.__dict__[key] = cf['para'].getboolean(key)
        else:
            try:
                pa.__dict__[key] = float(value)
            except:
                pa.__dict__[key] = value
    return pa

if __name__ == '__main__':
    cfg_file = '../Scripts/paraRF.cfg'
    parameter = readpara(cfg_file)
    print(parameter.get_para())
