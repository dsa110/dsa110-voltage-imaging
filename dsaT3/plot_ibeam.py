import sys

import glob
import numpy as np
import matplotlib.pylab as plt

import filplot_funcs

if __name__=='__main__':
#    /data/dsa110/T1/corr02/2021_9_1_18_33_28/fil_210902aabz/210902aabz_69.fil
    BASEDIR='/data/dsa110/T1/'
    datestr = sys.argv[1]
    candname = sys.argv[2]
    ibeam = int(sys.argv[3])
    ibox = int(sys.argv[4])
    dm = float(sys.argv[5])
    fpath = BASEDIR+'corr*/'+datestr+'/'+'fil_'+candname+'/*_'+str(ibeam)+'*.fil'
    fl = glob.glob(fpath)    
#    fn, dm, ibox = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
    fn = fl[0]
    fnout = fn.split('/')[-1]
    dataft, datadm, tsdm0, dms, datadm0 = filplot_funcs.proc_cand_fil(fn, dm, 1, snrheim=-1,
                                               pre_rebin=1, nfreq_plot=64,
                                               ndm=64, rficlean=True, norm=True)
    fnout = '/home/ubuntu/connor/data/%s_rficlean.npy' %  fn.split('/')[-1]
    print('Saving high res data to:\n%s' % fnout)
    np.save(fnout, dataft)
    dataft = dataft.downsample(ibox)
#    datadm = datadm.reshape(datadm.shape[0], -1, ibox).mean(-1)

    nfreq, ntime = dataft.shape
    xminplot,xmaxplot = 200,800 # milliseconds                                                                            
    tmin, tmax = 0., 1e3*dataft.header['tsamp']*ntime
    tarr = np.linspace(tmin, tmax, ntime)
    freqmax = dataft.header['fch1']
    freqmin = freqmax + dataft.header['nchans']*dataft.header['foff']

    dm_min, dm_max = dms.min(), dms.max()
    extentft=[tmin,tmax,freqmin,freqmax]
    extentdm=[tmin, tmax, dm_min, dm_max]
    
    fig = plt.figure(figsize=(8,10))
    plt.subplot(311)
    plt.imshow(dataft, aspect='auto', extent=extentft)
    plt.xlabel('Time (ms)')
    plt.ylabel('Freq ')
    plt.subplot(312)
    plt.plot(tarr, dataft.mean(0))
    plt.xlabel('Time (ms)')
    plt.xlim(0,tarr.max())
    plt.subplot(313)
    plt.imshow(datadm, aspect='auto', extent=extentdm)
    plt.suptitle('%s candname:%s \nDM:%0.1f boxcar:%d ibeam:%d' % (datestr, candname, dm, ibox, ibeam), color='C1')
    plt.show()

    
