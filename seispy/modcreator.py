import numpy as np
from scipy.interpolate import griddata, interp1d
from seispy.rfcorrect import DepModel
from seispy.setuplog import setuplog
import argparse


class ModCreator():
    def __init__(self, log=None):
        if log is None:
            self.logger = setuplog()
        else:
            self.logger = log
        self.lats = np.array([])
        self.lons = np.array([])
        self.deps = np.array([])
        self.vps = np.array([])
        self.vss = np.array([])
        self.ddvsddvp = np.array([[150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700],
                                 [1.85, 1.75, 1.75, 1.75, 1.79, 1.75, 2, 1.79, 1.61, 1.67, 1.67, 2.13]])

    def init_grid(self, latmin, latmax, lonmin, lonmax, depmin=0, depmax=850, horival=0.5, depval=50):
        self.logger.ModCreatorlog.info('Setting up grids.')
        self.lat_inter = np.arange(latmin, latmax, horival)
        self.lon_inter = np.arange(lonmin, lonmax, horival)
        self.dep_inter = np.arange(depmin, depmax, depval)
        self.dep_grid, self.lat_grid, self.lon_grid = np.meshgrid(self.dep_inter, self.lat_inter, self.lon_inter, indexing='ij')
        self.depmod = DepModel(self.dep_inter)

    def gridvel(self, **kwargs):
        self.logger.ModCreatorlog.info('Interpolating unstructured velocity to grid')
        points = np.array([self.deps, self.lats, self.lons]).T
        self.vp_grid = griddata(points, self.vps, (self.dep_grid, self.lat_grid, self.lon_grid), **kwargs)
        self.fill_1d('vp')
        if self.vss.size == 0:
            self.calc_vs()
        else:
            self.vs_grid = griddata(points, self.vss, (self.dep_grid, self.lat_grid, self.lon_grid), **kwargs)
            self.fill_1d('vs')

    def calc_vs(self):
        inter_ddvsddvp = interp1d(self.ddvsddvp[0], self.ddvsddvp[1], fill_value="extrapolate")(self.dep_inter)
        mesh_ddvsddvp, _, _ = np.meshgrid(inter_ddvsddvp, self.lat_inter, self.lon_inter, indexing='ij')
        vp_grid1d, _, _ = np.meshgrid(self.depmod.vp, self.lat_inter, self.lon_inter, indexing='ij')
        vs_grid1d, _, _ = np.meshgrid(self.depmod.vs, self.lat_inter, self.lon_inter, indexing='ij')
        self.dvp_grid = ((self.vp_grid - vp_grid1d) / vp_grid1d) * 100
        self.dvs_grid = self.dvp_grid * mesh_ddvsddvp
        self.vs_grid = vs_grid1d * (self.dvs_grid*0.01+1)
    
    def fill_1d(self, type='vp'):
        for i, _ in enumerate(self.dep_inter):
            self.__dict__[type+'_grid'][i][np.where(np.isnan(self.__dict__[type+'_grid'][i, :, :]))] = self.depmod.__dict__[type][i]

    def savenpz(self, filename):
        np.savez(filename, dep=self.dep_inter, lat=self.lat_inter, lon=self.lon_inter, vp=self.vp_grid, vs=self.vs_grid)

    @classmethod
    def read_txt(cls, datfile, usecols=(0, 1, 2, 3, 4), comments='#', delimiter=None):
        mc = cls()
        if len(usecols) == 5:
            mc.lats, mc.lons, mc.deps, mc.vps, mc.vss = np.loadtxt(datfile, unpack=True, usecols=usecols, comments=comments, delimiter=delimiter)
        elif len(usecols) == 4:
            mc.lats, mc.lons, mc.deps, mc.vps = np.loadtxt(datfile, unpack=True, usecols=usecols, comments=comments, delimiter=delimiter)
        else:
            raise ValueError('4 or 5 cols is valid in usecols')
        return mc


def veltxt2mod():
    parser = argparse.ArgumentParser(description="Convert velocity model with text format to npz format for 3-D moveout correction")
    parser.add_argument('velfile', help='velocity file in text format')
    parser.add_argument('-r', help='Range of the grid region with interval in degree',
                        type=str, required=True, metavar='lonmin/lonmax/latmin/latmax/val')
    parser.add_argument('-c', help='Which columns to read, with 0 being the first, order by lat/lon/dep/vp/vs. If 4 columns are accepted,'
                        'the fourth column represents Vp, and the Vs will be calculated as a relationship between'
                        'the Vp/Vs and the depth in Cammarano et al., (2003), defaults to 0/1/2/3/4',
                        type=str, default='0/1/2/3/4', metavar='ncol_lat/ncol_lon/ncol_dep/ncol_vp[/ncol_vs]')
    parser.add_argument('-d', help='Range of the grid region with interval in km, defaults to 0/850/50', 
                        type=str, default='0/850/50', metavar='depmin/depmax/val')
    parser.add_argument('-m', help='Method of interpolation. One of \'linear\', \'nearest\', \'cubic\', defaults to \'linear\'',
                        type=str, default='linear', metavar='linear|nearest|cubic')
    parser.add_argument('-o', help='output filename with \'npz\' format, defaults to ./mod3d.npz', default='./mod3d.npz')
    args = parser.parse_args()
    try:
        cols = [int(value) for value in args.c.split('/')]
    except:
        raise ValueError('Error in parsing -c option')
    try:
        gridrange = [float(value) for value in args.r.split('/')]
    except:
        raise ValueError('Error in parsing -r option')
    try:
        deprange = [float(value) for value in args.d.split('/')]
    except:
        raise ValueError('Error in parsing -d option')
    if args.m not in ['linear', 'nearest', 'cubic']:
        raise ValueError('The method of interpolation must be in \'linear\', \'nearest\', \'cubic\'')

    mc = ModCreator.read_txt(args.velfile, usecols=cols)
    mc.init_grid(gridrange[2], gridrange[3], gridrange[0], gridrange[1], depmin=deprange[0],
                 depmax=deprange[1], horival=gridrange[4], depval=deprange[2])
    mc.gridvel(method=args.m)
    mc.savenpz(args.o)


if __name__ == '__main__':
    veltxt2mod()

