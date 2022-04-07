"""Simple test of template ms writing."""

import os
import re
import subprocess
import datetime
from pkg_resources import resource_filename
import astropy.units as u
from dsacalib.ms_io import uvh5_to_ms
from dsavim.utils import load_params
from dsavim.voltages_to_ms import get_tstart_from_json
from dsavim.generate_uvh5 import generate_uvh5


PARAMFILE = resource_filename('dsaT3', 'data/T3_parameters.yaml')
# For continuum sources
NTINT = 512
NFINT = 8
START_OFFSET = 0
STOP_OFFSET = None

VOLTAGE_PATH = 'corr08_J1459+7140obw_data.out'
HEADER_PATH = 'J1459+7140obw.json'
DECLINATION_DEG = 71.67599222814322

class Logger():
    def __init__(self):
        current_time = datetime.datetime.now().strftime('%Y%m%dT%H%M%S')
        self.filepath = f'./test_template_writing_{current_time}_log.txt'
        with open(self.filepath, 'wt') as f:
            f.write(f'Starting at {current_time}\n')

    def info(self, message):
        with open(self.filepath, 'at') as f:
            f.write(f'\n{message}\n')

LOGGER = Logger()

def test_template_writing() -> None:
    """Generate template and data ms's using the same metadata.

    These are used in tests of writing the templates.
    """
    tstart = get_tstart_from_json(HEADER_PATH)
    correlation_parameters = get_correlation_parameters(NTINT)

    correlated_paths = dict({})
    LOGGER.info('Correlating data.')
    correlated_path = correlate_file(VOLTAGE_PATH, correlation_parameters)

    corrnode = re.findall('corr\d\d', correlated_path)[0]
    correlated_paths[corrnode] = correlated_path
    correlated_paths_string = '\n'.join([
        f'{key} {value}' for key, value in correlated_paths.items()])
    LOGGER.info(f'Correlated_files:\n{correlated_paths_string}')

    LOGGER.info('Making template files.')
    generate_ms('template', tstart, correlated_paths, correlation_parameters)

    LOGGER.info('Making data files.')
    generate_ms('data', tstart, correlated_paths, correlation_parameters)

    LOGGER.info('Finished.')

def generate_ms(candname: str, tstart: "astropy.time.Time", correlated_paths: dict, correlation_parameters: dict) -> None:
    """Generate a measurement set."""
    outfile_name = f'{correlation_parameters["corrdir"]}/{candname}'
    # if os.path.exists(f'{outfile_name}_{corrnode}.hdf5'):
    #     os.remove(f'{outfile_name}_{corrnode}.hdf5')

    uvh5_path = generate_uvh5(
        outfile_name,
        DECLINATION_DEG*u.deg,
        tstart,
        ntint=NTINT,
        nfint=NFINT,
        filelist=correlated_paths,
        start_offset=START_OFFSET,
        end_offset=STOP_OFFSET)
    LOGGER.info(f'Created {uvh5_path}.')

    # Generate ms
    outfile_name = f'{correlation_parameters["msdir"]}/{candname}'
    uvh5_to_ms([uvh5_path], outfile_name, fringestop=False)
    LOGGER.info(f'Created {outfile_name}.ms.')

def get_correlation_parameters(ntint: int) -> dict:
    input_params = load_params(PARAMFILE)
    deltat_ms = ntint*input_params['deltat_s']*1e3
    nants = len(input_params['antennas'])
    corrdir = input_params['corrdir']

    params = {
        'nants': nants,
        'deltat_ms': deltat_ms,
        'deltaf_MHz': input_params['deltaf_MHz'],
        'corrdir': corrdir,
        'msdir': input_params['msdir']
    }

    return params

def correlate_file(voltage_path: str, params: dict) -> str:
    """Correlate a voltage file for a single corr node."""
    output_path = f'{voltage_path}.corr'

    if not os.path.exists('{0}.corr'.format(voltage_path)):
        command = (
            '/home/ubuntu/proj/dsa110-shell/dsa110-bbproc/dsacorr '
            f'-d {voltage_path} -o {output_path} -t {params["deltat_ms"]} '
            f'-f {params["deltaf_MHz"]} -a {params["nants"]}')
        LOGGER.info(f'Running command:\n{command}')

        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            shell=True)

        proc_stdout = str(process.communicate()[0].strip())
        LOGGER.info(proc_stdout)

    return output_path

if __name__=='__main__':
    test_template_writing()
