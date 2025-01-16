from seispy.eq import rotateZNE
from obspy import read
from obspy.core.stream import Stream

st = read()
st.plot()

for tr in st:
    tr.data = tr.data * (-1)

print(st[0].data)

st.plot()
