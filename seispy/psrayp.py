import os
from obspy.taup import TauPyModel
import numpy as np
from seispy.geo import srad2skm
from scipy.interpolate import interpn
import time
import subprocess

def makepheaselist(layers):
    ph_list = ['P'+str(int(depth))+'s' for depth in layers]
    return ph_list


def PsRayp_test(layers, dist, dep):
    model = TauPyModel(model="iasp91")
    ph_list = makepheaselist(layers)
    t1 = time.time()
    arrs = model.get_ray_paths(dist, dep, ph_list)
    print(time.time()-t1)
    arr_num = len(arrs)
    rayp_list = np.zeros([arr_num, 2])
    for i in range(arr_num):
        rayp_list[i][0] = arrs[i].ray_param
        rayp_list[i][1] = int(arrs[i].name.strip('Ps'))
    rayp_list.sort(axis=0)
    print(time.time() - t1)
    return(rayp_list)


class PsRayp(object):
    def __init__(self, dis, dep, laymin=0, laymax=800):
        self.dis = dis
        self.dep = dep
        self.layers = np.arange(laymin, laymax)
        self.real_layers = np.array([])
        self.fake_layers = np.array([])
        self.real_idx = np.array([])
        self.fake_idx = np.array([])
        self.rayp = np.zeros((len(dis), len(dep), len(self.layers)))

    def make_phase_list(self):
        self.real_idx = np.where(self.layers >= 11)[0]
        self.fake_idx = np.where(self.layers < 11)[0]
        if self.real_idx.size == 0:
            raise ValueError('Max layer must greater than 8 km')
        self.real_layers = self.layers[self.real_idx]
        self.fake_layers = self.layers[self.fake_idx]
        with open('/tmp/tmp_phlst.txt', 'w+') as f:
            for lay in self.real_layers:
                f.write('P{}s\n'.format(lay))

    def taup_rayp(self, this_dis=50, this_dep=10, taup='taup_time'):
        """
        :param taup: Path to taup_time
        :return:
        """
        if not os.path.exists('/tmp/tmp_phlst.txt'):
            raise FileNotFoundError('Please excute \'make_phase_list\' first')
        cmd = '{} -h {} -pf {} -deg {} --rayp'.format(taup, this_dep, '/tmp/tmp_phlst.txt', this_dis)
        s = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        # s.communicate(cmd.encode())
        rayp = np.array(s.stdout.read().decode().strip().split())
        out_rayp = np.zeros([2, self.layers.shape[0]])
        for lay in self.fake_idx:
            out_rayp[0, lay] = self.layers[lay]
            out_rayp[1, lay] = rayp[0]
        for lay in self.real_idx:
            out_rayp[0, lay] = self.layers[lay]
            out_rayp[1, lay] = rayp[lay-11]
        return out_rayp

    def get_rayp(self):
        for i in range(self.dis.shape[0]):
            print('{}'.format(self.dis[i]))
            for j in range(self.dep.shape[0]):
                self.rayp[i, j, :] = self.taup_rayp(this_dis=self.dis[i], this_dep=self.dep[j])[1]

    def save(self, path='Ps_rayp'):
        np.savez(path, dis=self.dis, dep=self.dep, layers=self.layers, rayp=self.rayp)


def gen_rayp_lib(out_path, dis, dep, laymin=0, laymax=800):
    if not (isinstance(dis, np.ndarray) and isinstance(dep, np.ndarray)):
        raise TypeError('dis and dep must be np.ndarray type')
    pr = PsRayp(dis, dep, laymin=laymin, laymax=laymax)
    pr.make_phase_list()
    pr.get_rayp()
    pr.save(path=out_path)


def get_psrayp(rayp_lib, dis, dep, layers):
    x_layers = np.zeros([len(layers), 3])
    for i in range(len(layers)):
        x_layers[i] = np.array([dis, dep, layers[i]])
    return interpn((rayp_lib['dis'], rayp_lib['dep'], rayp_lib['layers']), rayp_lib['rayp'], x_layers, fill_value=None)


if __name__ == '__main__':
    layers = np.arange(11, 300)
    dep = np.arange(0, 300, 0.5)
    dis = np.arange(30, 90)
    rayp_path = '/Users/xumj/Researches/Ps_rayp.npz'
    rayp_lib = np.load(rayp_path)
    print(get_psrayp(rayp_lib, 50, 10, np.arange(0,100)))
