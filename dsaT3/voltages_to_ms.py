from dsaT3.utils import get_declination_mjd
from dsaT3.T3imaging import get_blen, generate_T3_ms, calibrate_T3ms

def __main__(name, filelist, ntint=96, start_offset=None, end_offset=None):
    # Get date from json file
    #
    # Get declination
    #
    # Correlate files
    #
    # Generate the measurement set
    msname = generate_T3_ms(
        name,
        pt_dec,
        tstart,
        ntint=ntint,
        nfint=1,
        filelist=filelist,
        start_offset=start_offset,
        end_offset=end_offset
    )
    print('{0} created'.format(msname))
