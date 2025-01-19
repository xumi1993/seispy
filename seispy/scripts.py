import numpy as np
import argparse
from seispy.rfcorrect import RFStation
from scipy.interpolate import interp1d
from seispy.utils import read_rfdep
from seispy.rf import RF
from seispy.recalrf import ReRF
import sys


def rfharmo():
    parser = argparse.ArgumentParser('Harmonic decomposition for extracting anisotropic and isotropic features from the radial and transverse RFs')
    parser.add_argument('rfpath', type=str, help="Path to PRFs")
    parser.add_argument('-e', help='Enlarge factor for plotting, defaults to 2', metavar='enf', default=2., type=float)
    parser.add_argument('-k', help="Suppress stacking RFs by back-azimuth, defaults to False", default=True,
                        action='store_false')
    parser.add_argument('-o', help="Specify output path for saving constant component.", metavar='outpath', default=None, type=str)
    parser.add_argument('-p', help='Figure output path, defaults to ./', metavar='figure_path', default='./', type=str)
    parser.add_argument('-s', help="Resample RFs with sampling interval of dt", metavar='dt', default=None, type=float)
    parser.add_argument('-t', help="Time window from tb to te for triming RFs, NOTE: do not insert space before this argument, defaults to -2/10",
                        metavar='tb/te', default='-2/10')
    args = parser.parse_args()
    rfsta = RFStation(args.rfpath)
    if args.s is not None:
        rfsta.resample(args.s)
    twin = [float(v) for v in args.t.split('/')]
    rfsta.harmonic(twin[0], twin[1])
    if args.o is not None:
        rfsta.harmo.write_constant(args.o)
    rfsta.harmo.plot(outpath=args.p, enf=args.e)


def rfani():
    parser = argparse.ArgumentParser(description="Estimate crustal anisotropy with a Joint inversion method. See Liu and Niu (2012) in detail.")
    parser.add_argument('rfpath', type=str, help="Path to PRFs")
    parser.add_argument('-t', help="Time window for searching Pms from tb to te", metavar='tb/te', required=True)
    parser.add_argument('-c', help="List file in text format for saving results, defaults to ./rfani.dat",
                        default="rfani.dat", metavar="list_file_name")
    parser.add_argument('-l', help="Half length of time window cut around Pms phase, defaults to 3s",
                        default=3, metavar="half_time_length", type=float)
    parser.add_argument('-r', help='Ray-parameter for moveout correction, defaults to 0.06 s/km',
                        default=0.06, metavar="rayp")
    parser.add_argument('-m', help='velocity model for moveout correction. \'iasp91\', \'prem\''
                        'and \'ak135\' is valid for internal model. Specify path to velocity model for the customized model.' 
                        'The format is the same as in Taup, but the depth should be monotonically increasing, defaults to \'iasp91\'',
                        default='iasp91', metavar="velocity_model")
    parser.add_argument('-o', dest='outpath', help="Directory to the image, defaults to current directory.", default='./')
    parser.add_argument('-p', help="If plot RFs stacked by back-azimuth, defaults to \'False\'",
                        dest="isplot", action='store_true', default=False)
    parser.add_argument('-w', help="Weights of 3 anisotropic methods (order by R cosine energy, R cross-correlation and T energy), defaults to 0.4/0.4/0.2",
                        dest='weight', default='0.4/0.4/0.2', metavar='w1/w2/w3')
    arg = parser.parse_args()
    weights = np.array(arg.weight.split('/')).astype(float)
    timewin = np.array(arg.t.split('/')).astype(float)
    rfsta = RFStation(arg.rfpath)
    bf, bt = rfsta.jointani(timewin[0], timewin[1], tlen=arg.l, weight=weights)
    if bf.size == 1:
        bf = np.append(bf, np.nan)
        bt = np.append(bt, np.nan)
    with open(arg.c, 'a+') as fid:
        fid.write('{}\t{:.3f}\t{:.3f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\n'.format(
                  rfsta.staname, rfsta.stla, rfsta.stlo, bf[0], bt[0], bf[1], bt[1]))
    if arg.isplot:
        rfsta.ani.plot_stack_baz(outpath=arg.outpath)
    rfsta.ani.plot_polar(outpath=arg.outpath)


def ccp3d():
    from seispy.ccp3d import CCP3D
    parser = argparse.ArgumentParser(description="3-D CCP stacking with spaced grid bins.")
    parser.add_argument('cfg_file', type=str, help='Path to CCP configure file')
    parser.add_argument('-s', help='Range for searching depth of D410 and D660, The results would be saved to \'peakfile\' in cfg_file',
                        metavar='d410min/d410max/d660min/d660max', default=None)
    arg = parser.parse_args()
    ccp = CCP3D(arg.cfg_file)
    ccp.initial_grid()
    ccp.stack()
    ccp.save_stack_data(ccp.cpara.stackfile)
    if arg.s:
        search_range = np.array(arg.s.split('/')).astype(float)
        ccp.search_good_410_660(*search_range)
        ccp.save_good_410_660(ccp.cpara.peakfile)


def ccp_profile():
    from seispy.ccpprofile import CCPProfile
    parser = argparse.ArgumentParser(description="Stack PRFS along a profile")
    parser.add_argument('cfg_file', type=str, help='Path to CCP configure file')
    parser.add_argument('-t', help='Output as a text file', dest='isdat', action='store_true')
    arg = parser.parse_args()
    if arg.isdat:
        typ = 'dat'
    else:
        typ = 'npz'
    ccp = CCPProfile(arg.cfg_file)
    ccp.initial_profile()
    ccp.stack()
    ccp.save_stack_data(format=typ)


def get_pierce_points():
    parser = argparse.ArgumentParser(description="Get pierce points with assumed depth")
    parser.add_argument('rfdepth_path', help="path to rfdepth file")
    parser.add_argument('-d', help="The depth in km", type=float, metavar='depth')
    parser.add_argument('-o', help="filename of output file, defaults to ./pierce_points.dat",
                        default='./pierce_points.dat', metavar='filename')
    arg = parser.parse_args()
    rfdep = read_rfdep(arg.rfdepth_path)
    if arg.d > rfdep[0]['depthrange'][-1]:
        parser.error('The depth exceed max depth in {}'.format(arg.rfdepth_path))
    with open(arg.o, 'w') as f:
        for sta in rfdep:
            for i in range(sta['piercelat'].shape[0]):
                la = interp1d(sta['depthrange'], sta['piercelat'][i])(arg.d)
                lo = interp1d(sta['depthrange'], sta['piercelon'][i])(arg.d)
                f.write('{:.4f} {:.4f}\n'.format(lo, la))
    

def common_parser():
    parser = argparse.ArgumentParser(description="Calculating RFs for single station")
    parser.add_argument('cfg_file', type=str, help='Path to RF configure file')
    parser.add_argument('-b', help='Correct back-azimuth. \nIf "baz" is specified, the corr_baz = raw_baz + baz. \n'
                                   'If there is no argument, the back-azimuth will be corrected with minimal '
                                   'energy of T component. The searching range is raw_baz +/- 90',
                                   dest='baz', nargs='?', const=0, type=float)
    parser.add_argument('-l', help="use local catalog, defaults to false", dest='islocal', action='store_true')
    parser.add_argument('-p', help='Wether or not manually pick arrival time and waveforms arround P/S phase with a GUI.',
                        action='store_true', default=False)
    parser.add_argument('-r', help='Reverse components: N, E or NE', dest='comp',
                        metavar='N|E|NE', default=None, type=str)
    parser.add_argument('-s', help='Switch the East and North components', dest='isswitch', action='store_true')
    parser.add_argument('-w', help='Write project to localfile', action='store_true')
    
    return parser


def parse_common_args(args):
    if args.comp is not None:
        args.comp = args.comp.upper()
        if args.comp == 'NE' or args.comp == 'EN':
            reverseE = True
            reverseN = True
        elif args.comp == 'E':
            reverseE = True
            reverseN = False
        elif args.comp == 'N':
            reverseE = False
            reverseN = True
        else:
            raise ValueError('component name must be in EN, E or N')
    else:
        reverseN = False
        reverseE = False
    return reverseE, reverseN


def prf():
    parser = common_parser()
    parser.add_argument('-f', help='Specify finallist for re-calculating RFs and -l is invalid in this pattern',
                        metavar='finallist', default=None)
    arg = parser.parse_args()
    if arg.f is not None:
        arg.islocal = False
        pjt = ReRF(arg.f, cfg_file=arg.cfg_file)
    else:
        pjt = RF(cfg_file=arg.cfg_file)
    if pjt.para.phase[0] != 'P':
        pjt.para.phase = 'P'
    pjt.para.switchEN = arg.isswitch
    pjt.para.reverseE ,pjt.para.reverseN= parse_common_args(arg)
    pjt.load_stainfo()
    if arg.f is None:
        pjt.search_eq(local=arg.islocal)
    pjt.match_eq()
    pjt.channel_correct()
    pjt.detrend()
    pjt.filter()
    pjt.cal_phase()
    if arg.f is None:
        pjt.drop_eq_snr()
    if arg.baz is not None and arg.baz != 0:
        pjt.baz_correct(correct_angle=arg.baz)
    elif arg.baz is not None and arg.baz == 0:
        pjt.baz_correct()
    else:
        pass
    if arg.w:
        pjt.savepjt()
    pjt.rotate()
    if arg.p:
        pjt.pick(prepick=False)
    pjt.trim()
    pjt.deconv()
    pjt.saverf()
    if arg.f is not None:
        pjt.write_list()


def srf():
    parser = common_parser()
    parser.add_argument('-i', help='Wether grid search incidence angle',
                        action='store_true')
    arg = parser.parse_args()
    pjt = RF(cfg_file=arg.cfg_file)
    if pjt.para.phase[0] != 'S':
        pjt.para.phase = 'S'
    pjt.para.switchEN = arg.isswitch
    pjt.para.reverseE ,pjt.para.reverseN= parse_common_args(arg)
    pjt.load_stainfo()
    pjt.search_eq(local=arg.islocal)
    pjt.match_eq()
    pjt.channel_correct()
    pjt.detrend()
    pjt.filter()
    pjt.cal_phase()
    pjt.drop_eq_snr()
    if arg.baz is not None and arg.baz != 0:
        pjt.baz_correct(correct_angle=arg.baz)
    elif arg.baz is not None and arg.baz == 0:
        pjt.baz_correct()
    else:
        pass
    pjt.rotate(search_inc=arg.i)
    if arg.p:
        pjt.pick()
    if arg.w:
        pjt.savepjt()
    pjt.trim()
    pjt.deconv()
    pjt.saverf()


def plot_rt():
    parser = argparse.ArgumentParser(description="Plot R(Q)&T components for P receiver functions (PRFs)")
    parser.add_argument('rfpath', help='Path to PRFs with a \'finallist.dat\' in it', type=str, default=None)
    parser.add_argument('-c', help='prime component to plot, defaults to \'R\'', 
                         default='R', type=str, metavar='[R|Q]')
    parser.add_argument('-e', help='Enlargement factor, defaults to 3', dest='enf', type=float, default=3, metavar='enf')
    parser.add_argument('-k', help='The key to sort PRFs, avialible for \'event\', \'evla\', \'evlo\', \'evdp\','
                                    '\'dis\', \'bazi\', \'rayp\', \'mag\', \'f0\', defaults to \'bazi\'', metavar='key',
                        default='bazi', type=str)
    parser.add_argument('-o', help='Output path without file name, defaults to current path', dest='output', default='./', type=str, metavar='outpath')
    parser.add_argument('-t', help='Specify figure format. f = \'.pdf\', g = \'.png\', defaults to \'g\'',
                        dest='format', default='g', type=str, metavar='f|g')
    parser.add_argument('-x', help='The max time scale in sec, defaults to 30s', default=30, type=float, metavar='max_time')
    arg = parser.parse_args()
    if arg.format not in ('f', 'g'):
        raise ValueError('Error: The format must be in \'f\' and \'g\'')
    rfsta = RFStation(arg.rfpath, prime_comp=arg.c)
    rfsta.plotrt(enf=arg.enf, out_path=arg.output, key=arg.k, outformat=arg.format, xlim=[-2, arg.x])


def plot_r():
    parser = argparse.ArgumentParser(description="Plot R&T receiver functions")
    parser.add_argument('rfpath', help='Path to PRFs with a \'finallist.dat\' in it', type=str)
    parser.add_argument('-c', help='prime component to plot, defaults to \'R\'', 
                         default='R', type=str, metavar='[R|Q|L|Z]')
    parser.add_argument('-e', help='Enlargement factor, defaults to 6', dest='enf', type=float, default=6, metavar='enf')
    parser.add_argument('-k', help='The key to sort PRFs, avialible for \'event\', \'evla\', \'evlo\', \'evdp\','
                                    '\'dis\', \'bazi\', \'rayp\', \'mag\', \'f0\', defaults to \'bazi\'', metavar='key',
                        default='bazi', type=str)
    parser.add_argument('-o', help='Output path without file name, defaults to current path', dest='output', default='./', type=str, metavar='outpath')
    parser.add_argument('-t', help='Specify figure format. f = \'.pdf\', g = \'.png\', defaults to \'g\'',
                    dest='format', default='g', type=str, metavar='f|g')
    parser.add_argument('-x', help='The max time scale in sec, defaults to 85s', default=85, type=float, metavar='max_time')

    arg = parser.parse_args()
    if arg.format not in ('f', 'g'):
        parser.error('Error: The format must be in \'f\' and \'g\'')
    rfsta = RFStation(arg.rfpath, prime_comp=arg.c)
    rfsta.plotr(arg.output, enf=arg.enf, key=arg.k, xlim=[-2, arg.x], outformat=arg.format)


def get_events():
    from obspy import UTCDateTime
    from seispy.io import Query
    parser = argparse.ArgumentParser(description="Get seismic events from IRIS Web-Service")
    parser.add_argument('-b', help='Start time, e.g., 20210101, 20210101020304',
                        metavar='datetime', default=None)
    parser.add_argument('-c', help='Catalog type',
                        metavar='GCMT|NEIC PDE|ISC', default=None)
    parser.add_argument('-d', help='Radial geographic constraints with center point and distance range',
                        metavar='<lat>/<lon>/<minradius>/<maxradius>', default=None)
    parser.add_argument('-e', help='End time, e.g., 20210101, 20210101020304',
                        metavar='datetime', default=None)
    parser.add_argument('-m', help='Magnitude range, optional for max magnitude',
                        metavar='<minmagnitude>[/<maxmagnitude>]', default=None)
    parser.add_argument('-r', help='Box range with min and max latitude and logitude, omited when specify \'-d\' ',
                        metavar='<lat>/<lon>/<minradius>/<maxradius>', default=None)
    parser.add_argument('-p', help='Focal depth, optional for max depth', 
                        metavar='<mindepth>[/<maxdepth>]', default=None)
    arg = parser.parse_args()
    args = {}
    if arg.c is not None:
        args['catalog'] = arg.c
    if arg.p is not None:
        try:
            values = [float(value) for value in arg.p.split('/')]
        except:
            raise ValueError('Error format with focal depth')
        if len(values) == 1:
            args['mindepth'] = values[0]
        elif len(values) == 2:
            args['mindepth'] = values[0]
            args['maxdepth'] = values[1]
        else:
            raise ValueError('Error format with focal depth')
    if arg.m is not None:
        try:
            values = [float(value) for value in arg.m.split('/')]
        except:
            raise ValueError('Error format with magnitude')
        if len(values) == 1:
            args['minmagnitude'] = values[0]
        elif len(values) == 2:
            args['minmagnitude'] = values[0]
            args['maxmagnitude'] = values[1]
        else:
            raise ValueError('Error format with focal depth')
    if arg.r is not None:
        try:
            values = [float(value) for value in arg.r.split('/')]
        except:
            raise ValueError('Error format with box range')
        if len(values) == 4:
            args['minlongitude'] = values[0]
            args['maxlongitude'] = values[1]
            args['minlatitude'] = values[2]
            args['maxlatitude'] = values[3]
        else:
            raise ValueError('Error format with box range')
    elif arg.d is not None:
        try:
            values = [float(value) for value in arg.d.split('/')]
        except:
            raise ValueError('Error format with Radial geographic constraints')
        if len(values) == 4:
            args['latitude'] = values[0]
            args['longitude'] = values[1]
            args['minradius'] = values[2]
            args['maxradius'] = values[3]
        else:
            raise ValueError('Error format with radial geographic constraints')
    else:
        pass
    if arg.b is not None:
        try:
            args['starttime'] = UTCDateTime(arg.b)
        except:
            raise ValueError('-b: Error format with time string')
    if arg.e is not None:
        try:
            args['endtime'] = UTCDateTime(arg.e)
        except:
            raise ValueError('-e: Error format with time string')
    if args == {}:
        parser.print_usage()
        sys.exit(1)
    query = Query()
    query.get_events(**args)
    for _, row in query.events.iterrows():
        print('{} {:.2f} {:.2f} {:.2f} {:.1f} {}'.format(
              row.date.isoformat(), row.evla, row.evlo,
              row.evdp, row.mag, row.magtype)
        )

def get_stations():
    from obspy import UTCDateTime
    from seispy.io import Query
    parser = argparse.ArgumentParser(description="Get stations from IRIS Web-Service")
    parser.add_argument('-S', help='Server name, defaults to IRIS', metavar='server', default='IRIS')
    parser.add_argument('-b', help='Start time, e.g., 20210101, 20210101020304', metavar='datetime', default=None)
    parser.add_argument('-c', help='Channel, wildcard is available like *,?,[EB]...', metavar='channel', default=None)
    parser.add_argument('-d', help='Radial geographic constraints with center point and distance range',
                        metavar='<lat>/<lon>/<minradius>/<maxradius>', default=None)
    parser.add_argument('-e', help='End time, e.g., 20210101, 20210101020304',
                        metavar='datetime', default=None)
    parser.add_argument('-n', help='Network name', metavar='network', default=None)
    parser.add_argument('-r', help='Box range with min and max latitude and logitude, omited when specify \'-d\'',
                        metavar='<lon1>/<lon2>/<lat1>/<lat2>', default=None)
    parser.add_argument('-s', help='Station name', metavar='station', default=None)
    arg = parser.parse_args()
    args = {}
    if arg.n is not None:
        args['network'] = arg.n
    if arg.s is not None:
        args['station'] = arg.s
    if arg.r is not None:
        try:
            values = [float(value) for value in arg.r.split('/')]
        except:
            raise ValueError('Error format with box range')
        if len(values) == 4:
            args['minlongitude'] = values[0]
            args['maxlongitude'] = values[1]
            args['minlatitude'] = values[2]
            args['maxlatitude'] = values[3]
        else:
            raise ValueError('Error format with box range')
    elif arg.d is not None:
        try:
            values = [float(value) for value in arg.d.split('/')]
        except:
            raise ValueError('Error format with Radial geographic constraints')
        args['latitude'] = values[0]
        args['longitude'] = values[1]
        args['minradius'] = values[2]
        args['maxradius'] = values[3]
    else:
        pass
    if arg.b is not None:
        try:
            args['starttime'] = UTCDateTime(arg.b)
        except:
            raise ValueError('-b: Error format with time string')
    if arg.e is not None:
        try:
            args['endtime'] = UTCDateTime(arg.e)
        except:
            raise ValueError('-e: Error format with time string')
    if arg.c is not None:
        args['channel'] = arg.c
    if args == {}:
        parser.print_usage()
        sys.exit(1)
    query = Query(arg.S)
    query.get_stations(**args)
    for net in query.stations:
        for sta in net:
            print('{} {} {:.6f} {:.6f} {:.4f} {} {}'.format(net.code, sta.code,
                  sta.latitude, sta.longitude, sta.elevation, sta.start_date, sta.end_date,
                  sta.restricted_status))