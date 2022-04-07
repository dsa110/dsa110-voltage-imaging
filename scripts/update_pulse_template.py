from pkg_resources import resource_filename
from dsavim.utils import load_params
from dsavim.update_template import update_template

paramfile = resource_filename("dsavim", "data/voltage_corr_parameters.yaml")
params = load_params(paramfile)

template_ms = f"{params['msdir']}/220121aaat.ms"
uvh5files = [f"{params['corrdir']}/220121aaat_{corr}.hdf5" for corr in params['ch0']]

update_template(template_ms, uvh5files)
