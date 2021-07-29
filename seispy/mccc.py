import numpy as np


def mccc(seis, dt, twin=0):
    data = np.zeros([len(seis), seis[0].data.shape[0]])
    for i in range(len(seis)):
        data[i, :] = seis[i].data
    nt = data.shape[1] * 2
    ns = data.shape[0]
    mask = np.ones(nt)
    tcc = np.zeros([ns, ns])
    tdel = np.zeros(ns)
    fft_all = np.fft.fft(data, nt, 1)
    fft_conj = fft_all.conj()

    if twin != 0:
        itw = np.fix(twin/(2*dt))
        mask = np.zeros(nt)
        mask[0:itw] = 1.0
        mask[nt-itw-1:nt] = 1.0
    for i in range(ns-1):
        ffis = np.conj(fft_conj[i])
#        ffis = np.conj(np.fft.fft(seis[i].data, nt))
        for j in np.arange(i+1, ns):
            ffjs = fft_all[j]
#            ffjs = np.fft.fft(seis[j].data, nt)
            ccf = np.fft.ifft(ffis*ffjs, nt).real*mask
            tcc[i, j] = np.argmax(ccf)

    (row, col) = np.where(tcc > nt/2)
    for i in range(row.shape[0]):
        tcc[row[i], col[i]] -= (nt+1)

    tcc *= dt

    for i in np.arange(0, ns):
        tdel[i] = (-np.sum(tcc[0:i+1, i]) + np.sum(tcc[i, i+1:ns]))/ns

    return tdel