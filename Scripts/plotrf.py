import os
import numpy as np
import re
try:
    import configparser
    config = configparser.ConfigParser()
except:
    import ConfigParser
    config = ConfigParser.ConfigParser()


def plot_R(rffiles, filenames, image_path, staname):
    fid_label = open("tmp_label", "w+")
    fid_baz = open("tmp_baz", "w+")
    fid_rf = open("tmp_rf", "w+")
    bazi = []
    evt_num = len(rffiles)
    gauss = rffiles[0].stats.sac.user1
    if gauss > 1.5:
        ax_len = 30
    else:
        ax_len = 80
    if not os.path.exists(image_path):
        os.makedirs(image_path)
    out_path = os.path.join(image_path, staname+'_R.ps')
    for i in range(evt_num):
        evtname = os.path.basename(filenames[i])
        evtname = re.split('[_|.]\w[_|.]',evtname)[0]
        timeaxis = rffiles[i].times()+rffiles[i].stats.sac.b
        bazi.append(rffiles[i].stats.sac.baz)
        fid_label.write('%d a %s\n' % (i+1, evtname))
        fid_baz.write('%d a %5.2f\n' % (i+1, bazi[i]))
        fid_rf.write('>\n')
        for j in range(rffiles[i].stats.npts):
            fid_rf.write('%6.3f %d %10.8f\n' % (timeaxis[j], i+1, rffiles[i].data[j]))
    fid_label.close()
    fid_baz.close()
    fid_rf.close()
    gmt = open('Plot_gmt.sh', 'w+')
    gmt.write('ps='+out_path+'\n')
    gmt.write('gmt gmtset FONT_ANNOT_PRIMARY 8p\n')
    gmt.write('gmt gmtset MAP_ANNOT_OFFSET_PRIMARY 3p\n')
    gmt.write('gmt gmtset MAP_GRID_PEN_PRIMARY 0.5p,gray\n')
    gmt.write('gmt gmtset MAP_TITLE_OFFSET 8p\n')
    gmt.write('gmt gmtset FONT_TITLE 13p\n')
    gmt.write('gmt gmtset FONT_LABEL 12p\n')
    gmt.write('gmt psbasemap -R-2/%d/0/%d -JX3i/6.5i -Bx5f1g5+l"Time after P (s)" -Byctmp_label -BWSen+t"R component (%s)" -K -X2.5i > $ps\n' % (ax_len, evt_num+2, staname))
    gmt.write('gmt pswiggle tmp_rf -W0.1p,gray -R -J -Z0.45 -G+red -G-blue -O -K >>$ps\n')
    gmt.write('gmt psxy -R -J -O -K -W1p >> $ps <<eof\n')
    gmt.write('0 0\n0 %d\neof\n' % (evt_num+2))
    gmt.write('gmt psbasemap -R0/360/0/%d -JX2i/6.5i -Bx60f10g60+l"Backazimuth (\\260)" -Byctmp_baz -BWSen -O -K -X3.7i>> $ps\n' % (evt_num+2))
    gmt.write("awk '{print $3,$1}' tmp_baz|gmt psxy -R -J -O -K -Sc0.05i -W0.1p,DODGERBLUE3 -GDODGERBLUE3 >> $ps\n")
    gmt.write('rm gmt*\n')
    gmt.write("rm tmp_label tmp_rf tmp_baz\n")
    gmt.close()
    os.system('bash Plot_gmt.sh')
    os.system("rm Plot_gmt.sh")

def plot_RT(rst, tst, filenames, image_path, staname):
    fid_label = open("tmp_label", "w+")
    fid_baz = open("tmp_baz", "w+")
    fid_rrf = open("tmp_rrf", "w+")
    fid_trf = open("tmp_trf", "w+")
    bazi = []
    evt_num = len(filenames)
    gauss = rst[0].stats.sac.user1
    if gauss > 1.5:
        ax_len = 30
    else:
        ax_len = 80
    if not os.path.exists(image_path):
        os.makedirs(image_path)
    out_path = os.path.join(image_path, staname+'_RT.ps')
    for i in range(evt_num):
        evtname = os.path.basename(filenames[i])
        evtname = re.split('[_|.]\w[_|.]',evtname)[0]
        timeaxis = rst[i].times()+rst[i].stats.sac.b
        bazi.append(rst[i].stats.sac.baz)
        fid_label.write('%d a %s\n' % (i+1, evtname))
        fid_baz.write('%d a %5.2f\n' % (i+1, bazi[i]))
        fid_rrf.write('>\n')
        fid_trf.write('>\n')
        for j in range(rst[i].stats.npts):
            fid_rrf.write('%6.3f %d %10.8f\n' % (timeaxis[j], i+1, rst[i].data[j]))
            fid_trf.write('%6.3f %d %10.8f\n' % (timeaxis[j], i+1, tst[i].data[j]))
    fid_label.close()
    fid_baz.close()
    fid_rrf.close()
    fid_trf.close()
    gmt = open('Plot_gmt.sh', 'w+')
    gmt.write('ps='+out_path+'\n')
    gmt.write('gmt gmtset FONT_ANNOT_PRIMARY 8p\n')
    gmt.write('gmt gmtset MAP_ANNOT_OFFSET_PRIMARY 3p\n')
    gmt.write('gmt gmtset MAP_GRID_PEN_PRIMARY 0.5p,gray\n')
    gmt.write('gmt gmtset MAP_TITLE_OFFSET 8p\n')
    gmt.write('gmt gmtset FONT_TITLE 13p\n')
    gmt.write('gmt gmtset FONT_LABEL 12p\n')
    gmt.write('gmt psbasemap -R-2/%d/0/%d -JX3i/6.5i -Bx5f1g5+l"Time after P (s)" -Byctmp_label -BWSen+t"R component (%s)" -K -X1.5i > $ps\n' % (ax_len, evt_num+2, staname))
    gmt.write('gmt pswiggle tmp_rrf -W0.1p,gray -R -J -Z0.45 -G+red -G-blue -O -K >>$ps\n')
    gmt.write('gmt psxy -R -J -O -K -W1p >> $ps <<eof\n')
    gmt.write('0 0\n0 %d\neof\n' % (evt_num+2))
    gmt.write('gmt psbasemap -R-2/%d/0/%d -JX3i/6.5i -Bx5f1g5+l"Time after P (s)" -BwSen+t"T component (%s)" -O -K -X3.2i >> $ps\n' % (ax_len, evt_num+2, staname))
    gmt.write('gmt pswiggle tmp_trf -W0.1p,gray -R -J -Z0.45 -G+red -G-blue -O -K >>$ps\n')
    gmt.write('gmt psxy -R -J -O -K -W1p >> $ps <<eof\n')
    gmt.write('0 0\n0 %d\neof\n' % (evt_num+2))
    gmt.write('gmt psbasemap -R0/360/0/%d -JX2i/6.5i -Bx60f10g60+l"Backazimuth (\\260)" -Byctmp_baz -BWSen -O -K -X3.5i>> $ps\n' % (evt_num+2))
    gmt.write("awk '{print $3,$1}' tmp_baz|gmt psxy -R -J -O -K -Sc0.05i -W0.1p,DODGERBLUE3 -GDODGERBLUE3 >> $ps\n")
    gmt.write('rm gmt*\n')
    gmt.write("rm tmp_label tmp_rrf tmp_trf tmp_baz\n")
    gmt.close()
    os.system('bash Plot_gmt.sh')
    os.system("rm Plot_gmt.sh")



