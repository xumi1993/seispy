"""
--------------------------------------------------------------------------------
Name:        get_cpt matplotlib colormap utility
Purpose:     an easy way to fetch .cpt colormap files, based on pycpt

Created:     2020.03
Copyright:   (c) Dimitrios Bouziotas (bouziot)
Licence:     GGNU General Public License v3 (GPL-3)
-You may freely copy, distribute and modify the software, in accordance with the provided license.
-You may not sublicense or hold the original author liable. This software comes with no warranty at all.
-Active contributions, forks and redevelopments are welcome.
-If you would like to include this software in your work, please reference it using the zenodo or github link. Please
also reference the original pycpt package (https://github.com/j08lue/pycpt)
--------------------------------------------------------------------------------
"""

__version__ = '0.1.0'
__copyright__ = """(c) 2020 Dimitrios Bouziotas"""
__license__ = "GGNU General Public License v3 (GPL-3)"

import os
import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

from urllib.request import urlopen

basedir = os.path.join(os.getcwd(),'cpt') # base path from where cpt files are fetched

def get_cmap(cpt_path, name=None, method='cdict', ret_cmap_type='LinearSegmented', N=None):
    """Get the cpt colormap as a LinearSegmented colormap. Utilizes the gmtColormap_openfile parser.
    Parameters
    ----------
    cpt_path : str, with one of the following options:
        - the full dir path to a .cpt file
        - the filename of a .cpt file that is in the local repo (check get_cpt.basedir)
        - a url.

    name : str, optional
        colormap name
        if not provided, the name will be derived automatically using _getname()

    method: str, optional
        Choose between 'cdict' and 'list'. The first one fetches all info from the .cpt file. The latter
        allows you to control the number of colors to sample from. Check gmtColormap_openfile for more info.

    N: int, optional
        the number of colors in the colormap to be returned. Can define the granularity of the returned colormap.
        Only useful when method='list'
    """
    # first get name
    if name is None:
        name = _getname(cpt_path)
    # case URL
    if 'http://' in cpt_path or 'https://' in cpt_path:
        with urlopen(cpt_path) as f:
            return gmtColormap_openfile(f, name=name, method=method, N=N, ret_cmap_type=ret_cmap_type)
    # case FULLPATH OR NAME
    else:
        with open(cpt_path) as f:
            return gmtColormap_openfile(f, name=name, method=method, N=N, ret_cmap_type=ret_cmap_type)

def get_listed_cmap(cpt_path, name=None, N=None):
    """Get the cpt colormap as a ListedColormap. Utilizes the gmtColormap_openfile parser.
    Parameters
    ----------
    cpt_path : str, with one of the following options:
        - the full dir path to a .cpt file
        - the filename of a .cpt file that is in the local repo (check get_cpt.basedir)
        - a url

    name : str, optional
        colormap name
        if not provided, the name will be derived automatically using _getname()

    N: int, optional
        the number of colors in the colormap to be returned. Leave None to derive the colors from the .cpt file.
        If you use a number less than the colors included in that file, a subset of colors will be returned.
    """
    # first get name
    if name is None:
        name = _getname(cpt_path)
    # case URL
    if 'http://' in cpt_path or 'https://' in cpt_path:
        with urlopen(cpt_path) as f:
            return gmtColormap_openfile(f, name=name, method='list', N=N, ret_cmap_type='Listed')
    # case FULLPATH OR NAME
    else:
        with open(cpt_path) as f:
            return gmtColormap_openfile(f, name=name, method='list', N=N, ret_cmap_type='Listed')

def _getname(cpt_path):
    """Internal function, fetches the name from a cpt filepath or url.
    Templates:
    'my.mby.cpt' -> 'my_mby' # NAME
    r'D:\matplotlib colormaps - cpt-city\cpt\mby.cpt' -> 'mby'  # FULLPATH
    'http://soliton.vm.bytemark.co.uk/pub/cpt-city/cmocean/haline.cpt'  -> 'haline' # URL
    """
    if 'http://' in cpt_path or 'https://' in cpt_path:  # CASE URL
        return '_'.join(cpt_path.split(r'/')[-1].split('.')[:-1])
    else:
        if os.path.exists(cpt_path):
            return cpt_path
        else:
            raise FileNotFoundError("CPT file {} not found".format(cpt_path))

def gmtColormap_openfile(cptf, name=None, method='cdict', N=None, ret_cmap_type='LinearSegmented'):
    """Read a GMT color map from an OPEN cpt file
    Edited by: bouziot, 2020.03

    Parameters
    ----------
    cptf : str, open file or url handle
        path to .cpt file

    name : str, optional
        name for color map
        if not provided, the file name will be used

    method : str, suggests the method to use.
    If method = 'cdict', generates the LinearSegmentedColormap using a color dictionary (cdict), disregarding any value in N.
    If method = 'list', generates the LinearSegmentedColor using the .from_list() method, passing a list of (value, (r,g,b)) tuples obtained from the cpt file. This allows the selection of colormap resolution by the user, using the N parameter

    N : int, the number of colors in the colormap. Only useful when method='list'.

    ret_cmap_type: str, the type of matplotlib cmap object to be returned. Accepts either 'LinearSegmented', which returns a matplotlib.colors.LinearSegmentedColormap, or 'Listed', which returns a ListedColormap
    In case 'Listed' is selected, the method argument from the user is ignored and method is set to 'list' ('Linear' doesn't work with 'cdict').
    N is then passed to matplotlib.colors.ListedColormap().
    - If N is set to None: all colors of the cpt file will be returned as a list.
    - In case of a user-defined N, colors will be truncated or extended by repetition (see matplotlib.colors.ListedColormap).

    Returns
    -------
    a matplotlib colormap object (matplotlib.colors.LinearSegmentedColormap or matplotlib.colors.ListedColormap)

    Credits
    -------
    This function originally appears in pycpt, extensive edits from bouziot, 2020.03
    Original work in: https://github.com/j08lue/pycpt
    LOG OF EDITS (2020.03):
        - Fixed bug when parsing non-split '#' lines in .cpt files
        - Fixed bug - not identifying the colorModel '#' line correctly
        - Fixed binary comparison performance (introduced in python 3)
        - Added functionality to return ListedColormaps and cmaps with custom colors (method, ret_cmap_type args)
        - Added global name handling externally (_getname() func)
    """
    methodnames = ['cdict', 'list'] # accepted method arguments
    ret_cmap_types = ['LinearSegmented', 'Listed', 'raw'] # accepted return matplotlib colormap types

    # generate cmap name
    if name is None:
        name = _getname(cptf.name)
    #    name = '_'.join(os.path.basename(cptf.name).split('.')[:-1])

    # process file
    x = []
    r = []
    g = []
    b = []
    lastls = None
    for l in cptf.readlines():
        ls = l.split()

        # skip empty lines
        if not ls:
            continue

        # parse header info
        # this leads to mistakes in some files...
        # if ls[0] in ["#", b"#"]:  # '#' is not always separated from other letters in some cases...
        #    if ls[-1] in ["HSV", b"HSV"]:
        #        colorModel = "HSV"
        #    else:
        #        colorModel = "RGB"
        #    continue

        # byte comparison is not feasible in python 3
        if (isinstance(l, bytes) and l.decode('utf-8')[0] in ["#", b"#"]) or (isinstance(l, str) and l[0] in ["#", b"#"]):
            if ls[-1] in ["HSV", b"HSV"]:
                colorModel = "HSV"
                continue
            elif ls[-1] in ["RGB", b"RGB"]:
                colorModel = "RGB"
                continue
            else: # case rogue comment, ignore
                continue


        # skip BFN info
        if ls[0] in ["B", b"B", "F", b"F", "N", b"N"]:
            continue

        # parse color vectors
        x.append(float(ls[0]))
        r.append(float(ls[1]))
        g.append(float(ls[2]))
        b.append(float(ls[3]))

        # save last row
        lastls = ls

    # check if last endrow has the same color, if not, append
    if not ((float(lastls[5]) == r[-1]) and (float(lastls[6]) == g[-1]) and (float(lastls[7]) == b[-1])):
        x.append(float(lastls[4]))
        r.append(float(lastls[5]))
        g.append(float(lastls[6]))
        b.append(float(lastls[7]))

    x = np.array(x)
    r = np.array(r)
    g = np.array(g)
    b = np.array(b)

    if colorModel == "HSV":
        for i in range(r.shape[0]):
            # convert HSV to RGB
            rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360., g[i], b[i])
            r[i] = rr ; g[i] = gg ; b[i] = bb
    elif colorModel == "RGB":
        r /= 255.
        g /= 255.
        b /= 255.

    red = []
    blue = []
    green = []
    xNorm = (x - x[0])/(x[-1] - x[0])

    # return colormap
    if method == 'cdict' and ret_cmap_type == 'LinearSegmented':
        # generate cdict
        for i in range(len(x)):
            red.append([xNorm[i],r[i],r[i]])
            green.append([xNorm[i],g[i],g[i]])
            blue.append([xNorm[i],b[i],b[i]])
        cdict = dict(red=red,green=green,blue=blue)
        #return cdict
        return mcolors.LinearSegmentedColormap(name=name,segmentdata=cdict)

    elif method == 'list' and ret_cmap_type == 'LinearSegmented':
        # generate list of values in the form of (value, c)
        outlist = []
        for i in range(len(x)):
            tup = (xNorm[i], (r[i],g[i],b[i]))
            outlist.append(tup)

        if N and type(N) == int:
        #return outlist
            return mcolors.LinearSegmentedColormap.from_list(name, outlist, N=N)
        else:
            raise TypeError("Using the method 'list' requires you to set a number of colors N.")

    elif ret_cmap_type == 'Listed':
        # generate list of values and return it as ListedColormap
        # returns both colors and the normalized positions (pos) where colors change, in the form of two outputs pos, colors
        pos_out = []
        colors_out = []
        for i in range(len(x)):
            pos_out.append(xNorm[i]) # list of positions
            colors_out.append(mcolors.to_hex( (r[i],g[i],b[i]))) # list of colors
        # return pos, color pairs
        print(colors_out)
        if N and type(N) == int and N<=len(colors_out):
            pos_out = pos_out[:N] #truncate positions to N
            return pos_out, mcolors.ListedColormap(colors_out, name=name, N=N)
        elif N is None:
            return pos_out, mcolors.ListedColormap(colors_out, name=name)
        else:
            raise TypeError("N has to be a number of colors that is less than the actual colors found in the .cpt file (" + str(len(colors_out)) + " colors found).")
    elif ret_cmap_type == 'raw':
        pos_out = []
        colors_out = []
        for i in range(len(x)):
            pos_out.append(xNorm[i]) # list of positions
            colors_out.append((r[i]*255, g[i]*255, b[i]*255))
        return pos_out, colors_out
    else:
        raise TypeError("method has to be one of the arguments: " + str(methodnames) + " and ret_cmap_type has to be one of the arguments: " + str(ret_cmap_types))

def plot_cmaps(cmap_list, width=6, cmap_height=0.5, axes_off=False):
    """Plot a colormap or list of colormaps with their names.
    Parameters
    -------
    cmap_list (str, cmap object or list of cmap objects anr strings): a list of colormaps to plot, either as cmap objects OR as preinstalled matplotlib colormap strings
    width (float): width of plot
    cmap_height (float): height of each colormap in plot
    axes_off (bool): boolean to erase axes

    Returns
    -------
    a matplotlib figure object (matplotlib.figure.Figure)

    Credits
    -------
    This function originally appears in pycpt, slight edits from bouziot, 2020.03
    https://github.com/j08lue/pycpt
    http://matplotlib.org/examples/color/colormaps_reference.html
    """
    if not isinstance(cmap_list, list):
        cmap_list = [cmap_list] # make subscriptable

    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))

    fig, axes = plt.subplots(nrows=len(cmap_list), figsize=(width,cmap_height*len(cmap_list)))
    fig.subplots_adjust(top=1, bottom=0, left=0, right=0.9)

    if len(cmap_list) == 1:
        cmap = cmap_list[0]
        # CASE ONE CMAP
        if isinstance(cmap, str):
                cmap = plt.get_cmap(cmap)
        axes.imshow(gradient, aspect='auto', cmap=cmap)
        pos = list(axes.get_position().bounds)
        x_text = pos[0] + pos[2] + 0.02
        y_text = pos[1] + pos[3]/2.
        fig.text(x_text, y_text, cmap.name, va='center', ha='left', fontsize=12)
        if axes_off:
            axes.set_axis_off()
        return fig

    else:
        # CASE MULTIPLE CMAPS
        for i, cmap in enumerate(cmap_list):
            if isinstance(cmap, str):
                cmap = plt.get_cmap(cmap)
            axes[i].imshow(gradient, aspect='auto', cmap=cmap)
            pos = list(axes[i].get_position().bounds)
            x_text = pos[0] + pos[2] + 0.02
            y_text = pos[1] + pos[3]/2.
            fig.text(x_text, y_text, cmap.name, va='center', ha='left', fontsize=12)
        # Turn off *all* ticks & spines, not just the ones with colormaps.
        if axes_off:
            for ax in axes:
                ax.set_axis_off()
        return fig

if __name__ == '__main__':
    #tests
    # test 1: FULL PATH, LinearSegmented, method cdict
    cpt_path = r'D:\Users\bouzidi\Desktop\matplotlib colormaps - cpt-city\cpt'
    cpt_fullpath = os.path.join(cpt_path, 'mby.cpt')

    a = get_cmap(cpt_fullpath)
    print(a)

    # test 2: LOCAL FILE, CHANGE BASEDIR
    basedir = r'D:\Users\bouzidi\Desktop\matplotlib colormaps - cpt-city\test\new_ctp'
    print(basedir)
    myctp2 = 'purple-orange-d15.cpt'
    pos, b = get_listed_cmap(myctp2)

    # test 3: url
    myurl = 'http://soliton.vm.bytemark.co.uk/pub/cpt-city/km/purple-orange-d15.cpt'
    c = get_cmap(myurl)
    print(c)
    print(c.name)

    # test 4: plots
    fig = plot_cmaps([a,b,c])
    plt.show()