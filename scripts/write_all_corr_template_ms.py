"""A script for creating a multi-spw measurement set from a single-spw template."""

from shutil import copytree
import argparse
import re
from pkg_resources import resource_filename
import numpy as np
from casacore.tables import table
from casatasks import virtualconcat
from dsaT3.utils import load_params
from dsaT3.update_template import update_template

PARAMFILE = resource_filename('dsaT3', 'data/T3_parameters.yaml')

# TODO: update_template_ms was designed for an entire ms with many frequencies, but
# TODO: this is too difficult to copy, read in etc.  Instead we will update each submms
# TODO: and then virtual concat them.  Calibration in this mode needs to be tested,
# TODO: and set_spectral_window should be moved to update_template_ms

def write_template(args: argparse.Namespace) -> None:
    """Create a template multi-spectral window ms."""
    params = get_parameters(args)

    template_ms = f'{params["msdir"]}/{params["template_name"]}.ms'

    for corr_node, uvh5file in params["uvh5files"].items():
        corr_freq = params["corr_node_freqs"][corr_node]
        corr_ms = f'{params["msdir"]}/{params["candname"]}_{corr_node}.ms'

        copytree(template_ms, corr_ms)
        set_spectral_window(corr_ms, corr_freq)
        update_template(corr_ms, [uvh5file])

    single_corr_ms_list = [f'{params["msdir"]}/{params["candname"]}_{corr_node}.ms'
                           for corr_node in params["uvh5files"]]
    multi_corr_ms = f'{params["msdir"]}/{params["candname"]}.ms'
    virtualconcat(single_corr_ms_list, multi_corr_ms)

def get_parameters(args: argparse.Namespace) -> dict:
    """Parse and calculate parameters needed for creating the ms."""
    all_params = load_params(PARAMFILE)

    # Set the frequencies for each corr node
    fobs = all_params['f0_GHz']+all_params['deltaf_MHz']*1e-3*np.arange(all_params['nchan'])
    corr_node_freqs = {}
    for corr, ch0 in all_params['ch0'].items():
        corr_node_freqs[corr] = fobs[ch0:(ch0+all_params['nchan_corr'])]

    # Set the filename for each corr node
    if args.uvh5file_dir:
        uvh5file_dir = args.uvh5file_dir
    else:
        uvh5file_dir = all_params['corrdir']

    uvh5files = {}
    for uvh5file in args.uvh5files:
        corr_node = re.findall('corr\d\d', uvh5file)[0]
        uvh5files[corr_node] = f'{uvh5file_dir}/{uvh5file}'

    params = {
        'msdir': all_params['msdir'],
        'template_name': args.template_name,
        'candname': args.candname,
        'corr_node_freqs': corr_node_freqs,
        'uvh5files': uvh5files
    }

    return params

def set_spectral_window(corr_ms: str, corr_freq: np.ndarray) -> None:
    """Set the spectral window parameters for the new ms.

    If the template spectral window has fewer channels than the corr node,
    a median is used to determine the output channels on the coarser reoslution.
    """
    with table(f'{corr_ms}/SPECTRAL_WINDOW', readonly=False) as tb:
        nchan_ms = np.array(tb.NUM_CHAN[:])[0]
        freq_out = np.median(corr_freq.reshape(nchan_ms, -1), axis=-1)
        tb.putcol('CHAN_FREQ', freq_out.reshape(1, nchan_ms))
        tb.putcol('REF_FREQUENCY', np.array([freq_out[0]]))

def parse_command_line_arguments() -> argparse.Namespace:
    """Parse commandline arguments."""
    parser = argparse.ArgumentParser(
        description='Make a multi-spectral window template ms.')
    parser.add_argument(
        '--template_name',
        type=str,
        help='name of the single spectral window ms to use a template',
        default='template_singlecorr')
    parser.add_argument(
        '--candname',
        type=str,
        help='name of the candidate for which to create the ms')
    parser.add_argument(
        '--uvh5files',
        type=str,
        nargs='+',
        help='the uvh5files to copy visibilities from')
    parser.add_argument(
        '--uvh5file_dir',
        type=str,
        help='the directory of the uvh5files',
        default='')

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    ARGS = parse_command_line_arguments()
    write_template(ARGS)
