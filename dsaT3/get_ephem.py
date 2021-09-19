import sys

import utils

if __name__=='__main__':
    if len(sys.argv)==1:
        ra,dec,mjd = utils.get_pointing_now()
    elif len(sys.argv)==2:
        mjd = float(sys.argv[1])
        ra,dec = utils.get_pointing_mjd(mjd)


    utils.get_galcoord(ra.value, dec.value)
    
    ind_near = utils.match_pulsar(ra, dec, thresh_deg=3.5)

    print("\n\nMJD: %0.5f"%mjd)
    print("RA and Dec: %0.2f %0.2f" % (ra.value,dec.value))
    if len(ind_near)==0:
        print("There are no pulsars within 3.5deg of beam center")
    else:
        print("There is/are %d pulsar(s) within 3.5deg of beam center:"%len(ind_near))

    for ii in ind_near:
        name_psrb = utils.query['PSRB'][ii]
        name_psrj = utils.query['PSRJ'][ii]
        dm = utils.query['DM'][ii]
        flux = utils.query['S1400'][ii]
        print('    %s/%s with DM=%0.1f pc cm**-3 and S1400=%0.3fmJy'%(name_psrb, name_psrj, dm, flux))
    

