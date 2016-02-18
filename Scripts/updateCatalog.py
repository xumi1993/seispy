#!/usr/bin/env python


try:
    import urllib.request as rq
except:
    import urllib as rq
import re
import sys, getopt, os
import datetime

def Usage():
    print('Usage: python updateCatalog.py -I[oldcatalog] -O[newcatalog]')
    print('-I   -- Old catalog')
    print('-O   -- New catalog')

try:
    opts,args = getopt.getopt(sys.argv[1:], "I:O:")
except:
    print('Arguments are not found!')
    Usage()
    sys.exit(1)
if opts == []:
    Usage()
    sys.exit(1)

outlog = ''
for op, value in opts:
    if op == "-I":
        inlog = value
    elif op == "-O":
        outlog = value
if outlog == '':
    outlog = inlog

url = 'http://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_QUICK/qcmt.ndk'
try:
    print('Connecting to http://www.ldeo.columbia.edu/.../qcmt.ndk')
    response = rq.urlopen(url)
except:
    print('Could not connect to http://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_QUICK/qcmt.ndk')

html = str(response.read())
version = sys.version_info.major
if version == 3:
    find_re = re.compile(r'[A-Z]+\s\d+/\d+/\d+\s.+?\\n',re.DOTALL)
else:
    find_re = re.compile(r'[A-Z]+\s\d+/\d+/\d+\s.+?[\n]',re.DOTALL)
fid_old = open(inlog,'r')
all_old_log = fid_old.readlines()
old_log = all_old_log[-1]
old_time_end = datetime.datetime(int(old_log.split()[0]),int(old_log.split()[1]),int(old_log.split()[2]),int(old_log.split()[4]),int(old_log.split()[5]),int(old_log.split()[6]))
fid_old.close()

if inlog != outlog:
    fid_new = open(outlog,'w+')
    fid_new.writelines(all_old_log)
else:
    fid_new = open(outlog,'a+')
i=0
maxcol = int(os.popen('stty size').read().strip().split()[1])
maxsharp = int(maxcol*0.8)
print('Writing station info to '+outlog)
for info in find_re.findall(html):
    num = len(find_re.findall(html))
    j = '#'*int((i/num)*maxsharp)
    sys.stdout.write(str(int((i/num)*100+1))+'% ||'+j+'->'+"\r")
    sys.stdout.flush()
    i+=1
    year = int(info.split()[1].split('/')[0])
    mon = int(info.split()[1].split('/')[1])
    day = int(info.split()[1].split('/')[2])
    hour = int(info.split()[2].split(':')[0])
    min = int(info.split()[2].split(':')[1])
    sec = float(info.split()[2].split(':')[2])
    sec = round(sec)
    if sec == 60:
        sec = 59
    lat = float(info.split()[3])
    lon = float(info.split()[4])
    dep = info.split()[5]
    mw = info.split()[7]
    evt_time = datetime.datetime(year,mon,day,hour,min,sec)
    if old_time_end < evt_time:
        fid_new.write('%d %d %d %s %d %d %d %6.2f %6.2f %s %s\n' % (year,mon,day,evt_time.strftime('%j'),hour,min,sec,lat,lon,dep,mw))
print('\nNew catalog file has writen to '+outlog+'.')


