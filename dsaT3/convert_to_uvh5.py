from dsacalib.ms_io import uvh5_to_ms
import glob

for name in ['J1459+7140obw']:
    msname = '/media/ubuntu/data/dsa110/imaging/{0}'.format(name)
    uvh5files = sorted(glob.glob(
        '/media/ubuntu/ssd/data/{0}_corr??.hdf5'.format(name)
    ))
    print(msname)
    print(len(uvh5files))
    print(uvh5files)
    uvh5_to_ms(
        uvh5files[:8],
        msname
    )
