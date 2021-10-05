import json, os, glob, sys
from dsautils import dsa_store
ds = dsa_store.DsaStore()
datestring = ds.get_dict('/cnf/datestring')

T3root = '/media/ubuntu/data/dsa110/T3/'

# copy from T3 directory 
os.system('mkdir -p '+T3root+datestring)
os.system('cp /home/ubuntu/data/T3/* '+T3root+datestring)

fls = glob.glob(T3root+datestring+'/*.json')

## let me know which ones will be archived
#for fl in fls:
#    f = open(fl)
#    de = json.load(f)
#    print('Key <save> not in {0}'.format(fl))
#    # Skip corr node json files without the save key if OoD archives twice
#    if de.get('save', False):
#        print('Will save voltages for ',de['trigname'])

for fl in fls:

    f = open(fl)
    de = json.load(f)
    #if de.get('save', False):
    if de['label']=='astrophysical':

        for corr in ['corr03', 'corr04', 'corr05', 'corr06', 'corr07', 'corr08', 'corr10', 'corr11', 'corr12', 'corr14', 'corr15', 'corr16', 'corr18', 'corr19', 'corr21', 'corr22']:

            if de[corr+'_header'] is True:

                outfile_h = T3root + datestring + '/'+corr+'_'+de['trigname']+'_header.json'
                
                if not os.path.exists(outfile_h):
                    print('copying header '+corr+' '+de['trigname'])
                    os.system('scp '+corr+'.sas.pvt:./data/'+de['trigname']+'_header.json '+outfile_h)
            
            if de[corr+'_data'] is True:

                outfile_d = T3root + datestring + '/'+corr+'_'+de['trigname']+'_data.out'

                if not os.path.exists(outfile_d):
                    print('copying data '+corr+' '+de['trigname'])
                    os.system('scp '+corr+'.sas.pvt:./data/'+de['trigname']+'_data.out '+outfile_d)
            
