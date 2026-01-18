#!/usr/bin/env python3
"""
Script to split a calibration table by SPW (spectral window) into separate files.
Each file is named as: YYYYMMDD_HHMMSS_XXMHz.bcal
"""

from casacore import tables
import numpy as np
import os
import shutil
import argparse
from datetime import datetime

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Split a calibration table by SPW into separate files')
parser.add_argument('input_caltable', type=str, 
                    help='Path to the input calibration table')
parser.add_argument('-o', '--output-dir', type=str, default=None,
                    help='Output directory (default: ./out_cal)')
args = parser.parse_args()

# Path to the input calibration table
input_caltable = args.input_caltable

# Output directory
if args.output_dir:
    output_dir = args.output_dir
else:
    # Default: use base name of input table
    base_name = os.path.basename(input_caltable.rstrip('/'))
    output_dir = f'./out_cal_{base_name}'

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

tb = tables.table(input_caltable)

# Get table description for creating new tables
table_desc = tb.getdesc()

# Get data
time_data = tb.getcol('TIME')
spw_data = tb.getcol('SPECTRAL_WINDOW_ID')
unique_spws = sorted(np.unique(spw_data))

# Get date/time from first TIME value
# CASA TIME is in seconds since MJD epoch (1858-11-17 00:00:00 UTC)
# MJD epoch in Unix timestamp: -3506716800.0 seconds
mjd_epoch_unix = -3506716800.0
casa_time = time_data[0]
unix_timestamp = casa_time + mjd_epoch_unix
dt = datetime.utcfromtimestamp(unix_timestamp)
date_str = dt.strftime('%Y%m%d')
time_str = dt.strftime('%H%M%S')

# Open SPECTRAL_WINDOW table to get frequencies
spw_table = tables.table(input_caltable + '/SPECTRAL_WINDOW')
chan_freq = spw_table.getcol('CHAN_FREQ')  # Shape: (nspw, nchannels)
spw_table.close()

# Process each SPW
for spw in unique_spws:
    print(f"Processing SPW {spw}")
    spw = int(spw)
    
    # Get start frequency for this SPW (first channel, in MHz)
    start_freq_mhz = chan_freq[spw, 0] / 1e6
    freq_int = int(np.round(start_freq_mhz))
    
    # Create output filename
    output_filename = f"{date_str}_{time_str}_{freq_int}MHz.bcal"
    output_path = os.path.join(output_dir, output_filename)
    
    # Remove existing output if it exists
    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    
    # Get rows for this SPW
    spw_mask = (spw_data == spw)
    spw_rows = np.where(spw_mask)[0]
    nrows_spw = len(spw_rows)
    
    if nrows_spw == 0:
        continue
    
    # Create new table with same structure
    tb_out = tables.table(output_path, table_desc, nrow=nrows_spw, ack=False)
    
    # Copy all columns for this SPW
    for colname in tb.colnames():
        try:
            col_data = tb.getcol(colname)
            # Select rows for this SPW
            if isinstance(col_data, np.ndarray):
                if len(col_data.shape) == 1:
                    col_data_spw = col_data[spw_mask]
                elif len(col_data.shape) == 2:
                    col_data_spw = col_data[spw_mask, :]
                elif len(col_data.shape) == 3:
                    col_data_spw = col_data[spw_mask, :, :]
                else:
                    col_data_spw = col_data[spw_mask, ...]
            else:
                # Handle list or other types
                col_data_spw = [col_data[i] for i in np.where(spw_mask)[0]]
                col_data_spw = np.array(col_data_spw)
            
            tb_out.putcol(colname, col_data_spw)
        except Exception:
            # Skip problematic columns
            continue
    
    # Set SPECTRAL_WINDOW_ID to 0 for all rows
    try:
        spw_id_data = tb_out.getcol('SPECTRAL_WINDOW_ID')
        spw_id_data[:] = np.int32(0)
        tb_out.putcol('SPECTRAL_WINDOW_ID', spw_id_data)
    except Exception:
        pass
    
    tb_out.flush()
    tb_out.close()
    
    # Copy subtables (ANTENNA, FIELD, OBSERVATION, SPECTRAL_WINDOW, HISTORY)
    subtable_names = ['ANTENNA', 'FIELD', 'OBSERVATION', 'SPECTRAL_WINDOW', 'HISTORY']
    
    for subtable_name in subtable_names:
        subtable_path = input_caltable + '/' + subtable_name
        output_subtable_path = output_path + '/' + subtable_name
        
        if os.path.exists(subtable_path):
            # For SPECTRAL_WINDOW, we only want the row for this SPW
            if subtable_name == 'SPECTRAL_WINDOW':
                # Copy the subtable
                if os.path.exists(output_subtable_path):
                    shutil.rmtree(output_subtable_path)
                os.makedirs(os.path.dirname(output_subtable_path), exist_ok=True)
                shutil.copytree(subtable_path, output_subtable_path)
                
                # Open and filter to only this SPW, then set it as SPW 0
                spw_tb = tables.table(output_subtable_path, readonly=False)
                if spw_tb.nrows() > spw:
                    # Keep only the row for this SPW (will become SPW 0)
                    # Read all data
                    spw_desc = spw_tb.getdesc()
                    spw_nrows = spw_tb.nrows()
                    
                    # Create a new table with only one row (SPW 0)
                    temp_path = output_subtable_path + '_temp'
                    if os.path.exists(temp_path):
                        shutil.rmtree(temp_path)
                    spw_tb_new = tables.table(temp_path, spw_desc, nrow=1, ack=False)
                    
                    # Copy the row for this SPW (from original SPW index)
                    for colname in spw_tb.colnames():
                        try:
                            col_data = spw_tb.getcol(colname)
                            if isinstance(col_data, np.ndarray):
                                if len(col_data.shape) == 1:
                                    spw_tb_new.putcol(colname, col_data[spw:spw+1])
                                elif len(col_data.shape) == 2:
                                    spw_tb_new.putcol(colname, col_data[spw:spw+1, :])
                                else:
                                    spw_tb_new.putcol(colname, col_data[spw:spw+1, ...])
                            else:
                                # Handle list or other types
                                spw_tb_new.putcol(colname, [col_data[spw]])
                        except Exception:
                            continue
                    
                    spw_tb_new.flush()
                    spw_tb_new.close()
                    spw_tb.close()
                    
                    # Replace the old subtable with the new one
                    shutil.rmtree(output_subtable_path)
                    shutil.move(temp_path, output_subtable_path)
                else:
                    spw_tb.close()
            else:
                # For other subtables, copy as-is
                if os.path.exists(output_subtable_path):
                    shutil.rmtree(output_subtable_path)
                os.makedirs(os.path.dirname(output_subtable_path), exist_ok=True)
                shutil.copytree(subtable_path, output_subtable_path)
    
    # Now update keywords with correct subtable paths (after subtables are copied)
    tb_out = tables.table(output_path, readonly=False)
    try:
        keywords = tb.getkeywords()
        if keywords:
            # Get absolute path for output
            abs_output_path = os.path.abspath(output_path)
            # Update subtable references to point to new table paths
            updated_keywords = {}
            for key, value in keywords.items():
                if isinstance(value, str) and 'Table:' in value:
                    # Update subtable path to point to new table (use absolute path)
                    subtable_name = key
                    updated_keywords[key] = f'Table: {abs_output_path}/{subtable_name}'
                else:
                    updated_keywords[key] = value
            tb_out.putkeywords(updated_keywords)
            tb_out.flush()
    except Exception as e:
        print(f"Warning: Could not update keywords: {e}")
    finally:
        tb_out.close()
    
    # Write table.info with Type and SubType
    table_info_path = os.path.join(output_path, 'table.info')
    with open(table_info_path, 'w') as f:
        f.write('Type = Calibration\n')
        f.write('SubType = B Jones\n')

# Close input table
tb.close()

print(f"Created {len(unique_spws)} files in {output_dir}")
