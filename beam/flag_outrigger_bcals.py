from casatools import table
from ovrolwasolar import flagging
import argparse
import os
import glob

def flag_outrigger(dataset, ref_ms):
    '''
    dataset: can be a MS or caltable
    ref_MS: MS from which the antenna names will be read
    
    Note that this function assumes that there is a single timeslice
    in the dataset. If there are more than one timelsice, this
    function will give wrong results.
    '''
    tb=table()
    core_ant_ids, exp_ant_ids = flagging.get_antids(ref_ms)
    tb.open(dataset,nomodify=False)
    try:
        flag=tb.getcol('FLAG')
        flag[:,:,exp_ant_ids]=True
        tb.putcol('FLAG',flag)
        tb.flush()
    finally:
        tb.close()
    return


# Parse command-line arguments
parser = argparse.ArgumentParser(description='Flag outrigger antennas in bcal files')
parser.add_argument('input_dir', type=str,
                    help='Input directory containing bcal files')
parser.add_argument('-ms', '--ref_ms', type=str,
                    default='/data07/peijinz/tests/Gregg_caltable/ms/20251124_192303_27MHz.ms',
                    help='Reference MS from which the antenna names will be read')
args = parser.parse_args()

input_dir = args.input_dir
beam_caltable_dir = f'{input_dir}_beam'

# Create output directory
os.makedirs(beam_caltable_dir, exist_ok=True)

# Find all bcal files in input directory
bcaltbs = glob.glob(os.path.join(input_dir, '*.bcal'))
bcaltbs.sort()

print(f"Found {len(bcaltbs)} bcal files in {input_dir}")
print(f"Output directory: {beam_caltable_dir}")

# Process each bcal file
for bcaltb in bcaltbs:
    bcaltb_bm = os.path.join(beam_caltable_dir, os.path.basename(bcaltb))
    print(f"Processing: {os.path.basename(bcaltb)}")
    
    # Copy bcal to beam directory
    os.system('cp -r ' + bcaltb + ' ' + bcaltb_bm)
    
    # Flag outriggers
    flag_outrigger(bcaltb_bm, args.ref_ms)
    print(f"  Flagged outriggers in: {bcaltb_bm}")

print(f"\nDone! Processed {len(bcaltbs)} bcal files")


