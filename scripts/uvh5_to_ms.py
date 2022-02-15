from dsacalib.ms_io import uvh5_to_ms
import glob

for name in ['220128aaas']:
    uvh5files = sorted(glob.glob(
        '/media/ubuntu/ssd/data/{0}_corr??.hdf5'.format(name)
    ))
    for uvh5file in uvh5files:
        msname = '/media/ubuntu/data/dsa110/imaging/{0}'.format(
            uvh5file.replace('.hdf5', '').split('/')[-1])
        uvh5_to_ms(
            uvh5file,
            msname
        )
