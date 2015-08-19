#!/usr/bin/env python

def Usage():
   print('./PlotR.py -Sstation -Iinpath -Ooutpath -Nnpts -Atimeaxis -Ttime_before/time_after')

import os
import numpy as np
import getopt
import sys
try:
    import configparser
    config = configparser.ConfigParser()
except:
    import ConfigParser
    config = ConfigParser.ConfigParser()

try:
   opts, args = getopt.getopt(sys.argv[1:], "S:I:O:N:T:A:")
except:
    print('Arguments are not found!')
    Usage()
    sys.exit(1)
if opts == []:
    Usage()
    sys.exit(1)

for op, value in opts:
    if op == "-S":
        staname = value
    elif op == "-I":
        in_path = value
    elif op == "-O":
        out_path = value
    elif op == "-N":
        RFlength = int(value)
    elif op == "-T":
        time_str = value
    elif op == "-A":
        axis_len = value
    else:
       Usage()
       sys.exit(1)

time_before = float(time_str.split('/')[0])
time_after = float(time_str.split('/')[1])
fid_list = open(os.path.join(in_path, staname+'finallist.dat'), 'r')
date_name = np.array([])
bazi = np.array([])
for line in fid_list.readlines():
    line = line.strip('\n')
    date_name = np.append(date_name, line.split()[0])
    bazi = np.append(bazi, float(line.split()[6]))

idx = np.argsort(bazi)
bazi = bazi[idx]
date_name = date_name[idx]
timeaxis = np.linspace(-time_before, time_after, RFlength, endpoint=True)
fid_tmp = open('tmp.lst', 'w+')
fid_yticklabel = open('yticklabel.txt', 'w+')
for i in range(len(bazi)):
    print(bazi[i], date_name[i])
    RF = np.loadtxt(os.path.join(in_path, date_name[i]+'_P_R.dat'))
    fid_yticklabel.write('%d a %s\n' % (i+1, date_name[i]))
    fid_tmp.write('>\n')
    for j in range(RFlength):
        fid_tmp.write('%6.3f %d %10.8f\n' % (timeaxis[j], i+1, RF[j]))
gmt = open('Plot_gmt.sh', 'w+')
gmt.write('ps='+out_path+'\n')
gmt.write('gmt psbasemap -R-2/%s/0/%d -JX4i/9i -Bx5f1+l"Time after P (s)" -Bycyticklabel.txt -BWSen+t"%s" -K -P -X2i > $ps\n' % (axis_len, len(bazi)+2, staname))
gmt.write('gmt pswiggle tmp.lst -W0.1p,gray -R -J -Z0.45 -G+red -G-blue -O -K >>$ps\n')
gmt.write('gmt psxy -R -J -O -K -W1p >> $ps <<eof\n')
gmt.write('0 0\n0 %d\neof\n' % (len(bazi)+2))
gmt.write('rm gmt* yticklabel.txt tmp.lst')
gmt.close()
os.popen('bash Plot_gmt.sh')
#os.popen('rm gmt* yticklabel.txt tmp.lst')
# gmt.write('open '+out_path)
