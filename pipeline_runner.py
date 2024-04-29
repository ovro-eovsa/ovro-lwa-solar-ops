import sys
sys.path.append('/data07/msurajit/ovro-lwa-solar-ops/')
import solar_realtime_pipeline as srp
import logging
from astropy.time import Time

"""
    Main routine of running the realtime pipeline. Example call
        pdsh -w lwacalim[00-09] 'conda activate suncasa && python /opt/devel/bin.chen/ovro-lwa-solar/operations/solar_realtime_pipeline.py 2023-11-21T15:50'
    Sometimes afer killing the pipeline (with ctrl c), one need to remove the temporary files and kill all the processes before restarting.
        pdsh -w lwacalim[00-09] 'rm -rf /fast/bin.chen/realtime_pipeline/slow_working/*'
        pdsh -w lwacalim[00-09] 'rm -rf /fast/bin.chen/realtime_pipeline/slow_slfcaled/*'
        pdsh -w lwacalim[00-09] 'pkill -u bin.chen -f wsclean'
        pdsh -w lwacalim[00-09] 'pkill -u bin.chen -f python'
"""

import argparse,os

parser = argparse.ArgumentParser(description='Solar realtime pipeline')
parser.add_argument('prefix', type=str, help='Timestamp for the start time. Format YYYY-MM-DDTHH:MM')
parser.add_argument('--end_time', default='2030-01-01T00:00', help='End time in format YYYY-MM-DDTHH:MM')
parser.add_argument('--interval', default=600., help='Time interval in seconds')
parser.add_argument('--nodes', default='0123456789', help='List of nodes to use')
parser.add_argument('--delay', default=60, help='Delay from current time in seconds')
parser.add_argument('--server', default=None, help='Name of the server where the raw data is located. Must be defined in ~/.ssh/config.')
parser.add_argument('--nolustre', default=False, help='If set, do NOT assume that the data are stored under /lustre/pipeline/ in the default tree', action='store_true')
parser.add_argument('--file_path', default='slow/', help='Specify where the raw data is located')
parser.add_argument('--proc_dir', default='/fast/bin.chen/realtime_pipeline/', help='Directory for processing')
parser.add_argument('--save_dir', default='/lustre/bin.chen/realtime_pipeline/', help='Directory for saving fits files')
parser.add_argument('--calib_dir', default='/lustre/bin.chen/realtime_pipeline/caltables/', help='Directory to calibration tables')
parser.add_argument('--calib_file', default='', help='Calibration file to be used yyyymmdd_hhmmss')
parser.add_argument('--alt_limit', default=15., help='Lowest solar altitude to start/end imaging')
parser.add_argument('--bmfit_sz', default=2, help='Beam fitting size to be passed to wsclean')
parser.add_argument('--do_refra', default=True, help='If True, do refraction correction', action='store_false')
parser.add_argument('--singlenode', default=False, help='If True, delay the start time by the node', action='store_true')
parser.add_argument('--logger_dir', default='/lustre/bin.chen/realtime_pipeline/logs/', help='Directory for logger files')
parser.add_argument('--logger_prefix', default='solar_realtime_pipeline', help='Prefix for logger file')
parser.add_argument('--logger_level', default=10, help='Specify logging level. Default to 10 (debug)')   
parser.add_argument('--keep_working_ms', default=False, help='If True, keep the working ms files after imaging', action='store_true')
parser.add_argument('--keep_working_fits', default=False, help='If True, keep the working fits files after imaging', action='store_true')  
parser.add_argument('--no_selfcal', default=False, help='If True, perform selfcal', action='store_true')
parser.add_argument('--no_imaging', default=False, help='If True, perform imaging', action='store_true')
parser.add_argument('--bands', '--item', action='store', dest='bands',
                    type=str, nargs='*', 
                    default=['32MHz', '36MHz', '41MHz', '46MHz', '50MHz', '55MHz', '59MHz', '64MHz', '69MHz', '73MHz', '78MHz', '82MHz'],
                    help="Examples: --bands 32MHz 46MHz 64MHz")

args = parser.parse_args()
if len(args.calib_file) == 15:
    calib_file = args.calib_file
else:
    print('Calibration tables not provided or recognized. Attempting to find those from default location on lwacalim.')
    calib_tables = glob.glob('/lustre/bin.chen/realtime_pipeline/caltables_latest/*.bcal')
    if len(calib_tables) > 10:
        calib_file = os.path.basename(calib_tables[0])[:15]
        print('Using calibration file {0:s}'.format(calib_file))
    else:
        print('No calibration files found. Abort.')
        sys.exit(0)

try:
    srp.run_pipeline(args.prefix, time_end=Time(args.end_time), time_interval=float(args.interval), nodes=args.nodes, \
                    delay_from_now=float(args.delay), server=args.server, lustre=(not args.nolustre), file_path=args.file_path,\
                    proc_dir=args.proc_dir, save_dir=args.save_dir, calib_dir=args.calib_dir, calib_file=calib_file, 
                    altitude_limit=float(args.alt_limit), logger_dir = args.logger_dir, logger_prefix=args.logger_prefix,\
                    logger_level=int(args.logger_level), do_refra=args.do_refra, 
                    multinode= (not args.singlenode), delete_working_ms=(not args.keep_working_ms), 
                    delete_working_fits=(not args.keep_working_fits), beam_fit_size=args.bmfit_sz, do_selfcal=(not args.no_selfcal),\
                    do_imaging=(not args.no_imaging), bands=args.bands)
except Exception as e:
    logging.error(e)
    raise e
