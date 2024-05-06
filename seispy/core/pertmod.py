from scipy.interpolate import interpn, interp1d
import numpy as np
from .depmodel import DepModel

class Mod3DPerturbation:
    """ 
    Perturbation model for 3D velocity model.
    """   
    def __init__(self, modpath, YAxisRange, velmod='iasp91'):
        dep_mod = DepModel(YAxisRange, velmod=velmod)
        self.model = np.load(modpath)
        new1dvp = interp1d(dep_mod.depths, dep_mod.vp)(self.model['dep'])
        new1dvs = interp1d(dep_mod.depths, dep_mod.vs)(self.model['dep'])
        new1dvp, _, _ = np.meshgrid(new1dvp, self.model['lat'], self.model['lon'], indexing='ij')
        new1dvs, _, _ = np.meshgrid(new1dvs, self.model['lat'], self.model['lon'], indexing='ij')
        self.dvp = (self.model['vp'] - new1dvp) / new1dvp
        self.dvs = (self.model['vs'] - new1dvs) / new1dvs
        self.cvp = dep_mod.vp
        self.cvs = dep_mod.vs

    def interpdvp(self, points):
        dvp = interpn((self.model['dep'], self.model['lat'], self.model['lon']), self.dvp, points,
                      bounds_error=False, fill_value=None)
        return dvp

    def interpdvs(self, points):
        dvs = interpn((self.model['dep'], self.model['lat'], self.model['lon']), self.dvs, points,
                      bounds_error=False, fill_value=None)
        return dvs
