from dsacalib.ms_io import uvh5_to_ms
import glob

for name in ["220128aaas"]:
    uvh5files = sorted(glob.glob(
        f"/media/ubuntu/ssd/data/{name}_corr??.hdf5"))
    for uvh5file in uvh5files:
        msname = f"/media/ubuntu/data/dsa110/imaging/{uvh5file.replace('.hdf5', '').split('/')[-1]}"
        uvh5_to_ms(uvh5file, msname)
