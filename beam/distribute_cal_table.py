#!/usr/bin/env python3
"""
Script to distribute calibration tables to various locations and remote nodes.
"""

import argparse
import os
import shutil
import glob
import subprocess

# Default paths
DEFAULT_SRC_CAL = "/data07/peijinz/tests/Gregg_caltable/out_cal_2025-12-05"
DEFAULT_SRC_CAL_BEAM = "/data07/peijinz/tests/Gregg_caltable/out_cal_2025-12-05_beam"

# Destination paths
DEST_CALTABLES = "/opt/devel/solarpipe/operation/caltab/caltables"
DEST_CALTABLES_LUSTRE = "/lustre/solarpipe/realtime_pipeline/caltables"
DEST_CALTABLES_BEAM = "/opt/devel/solarpipe/operation/caltab/caltables_beam"
DEST_CALTABLES_BEAM_LUSTRE = "/lustre/solarpipe/realtime_pipeline/caltables_beam"

DEST_CALTABLES_LATEST = "/opt/devel/solarpipe/operation/caltab/caltables_latest"
DEST_CALTABLES_BEAM_LATEST = "/opt/devel/solarpipe/operation/caltab/caltables_beam_latest"
DEST_CALTABLES_LATEST_LUSTRE = "/lustre/solarpipe/realtime_pipeline/caltables_latest"
DEST_CALTABLES_BEAM_LATEST_LUSTRE = "/lustre/solarpipe/realtime_pipeline/caltables_beam_latest"

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Distribute calibration tables to various locations')
parser.add_argument('--src-cal', type=str, default=DEFAULT_SRC_CAL,
                    help=f'Source directory for calibration tables (default: {DEFAULT_SRC_CAL})')
parser.add_argument('--src-cal-beam', type=str, default=DEFAULT_SRC_CAL_BEAM,
                    help=f'Source directory for beam calibration tables (default: {DEFAULT_SRC_CAL_BEAM})')
parser.add_argument('--skip-pdsh', action='store_true',
                    help='Skip pdsh commands for remote distribution')
args = parser.parse_args()

src_cal = args.src_cal
src_cal_beam = args.src_cal_beam

def copy_files(src_dir, dest_dir, description):
    """Copy all files from src_dir to dest_dir."""
    if not os.path.exists(src_dir):
        print(f"Warning: Source directory does not exist: {src_dir}")
        return
    
    os.makedirs(dest_dir, exist_ok=True)
    
    # Find all files/directories in source
    items = glob.glob(os.path.join(src_dir, '*'))
    
    if not items:
        print(f"Warning: No files found in {src_dir}")
        return
    
    print(f"\n{description}")
    print(f"  Copying from: {src_dir}")
    print(f"  Copying to: {dest_dir}")
    
    for item in items:
        dest_item = os.path.join(dest_dir, os.path.basename(item))
        try:
            if os.path.isdir(item):
                if os.path.exists(dest_item):
                    shutil.rmtree(dest_item)
                shutil.copytree(item, dest_item)
            else:
                if os.path.exists(dest_item):
                    os.remove(dest_item)
                shutil.copy2(item, dest_item)
            print(f"    Copied: {os.path.basename(item)}")
        except Exception as e:
            print(f"    Error copying {os.path.basename(item)}: {e}")

def remove_files(dest_dir, description):
    """Remove all files from dest_dir."""
    if not os.path.exists(dest_dir):
        print(f"Warning: Destination directory does not exist: {dest_dir}")
        return
    
    print(f"\n{description}")
    print(f"  Removing files from: {dest_dir}")
    
    items = glob.glob(os.path.join(dest_dir, '*'))
    for item in items:
        try:
            if os.path.isdir(item):
                shutil.rmtree(item)
            else:
                os.remove(item)
            print(f"    Removed: {os.path.basename(item)}")
        except Exception as e:
            print(f"    Error removing {os.path.basename(item)}: {e}")

def run_pdsh_command(cmd, description):
    """Run a pdsh command."""
    print(f"\n{description}")
    print(f"  Running: {cmd}")
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"    Success")
            if result.stdout:
                print(f"    Output: {result.stdout}")
        else:
            print(f"    Error (return code {result.returncode})")
            if result.stderr:
                print(f"    Error output: {result.stderr}")
    except Exception as e:
        print(f"    Exception: {e}")

print("="*60)
print("Distributing Calibration Tables")
print("="*60)

# Step 1: Copy src_cal/* to caltables directories
copy_files(src_cal, DEST_CALTABLES, "Step 1a: Copying calibration tables to caltab/caltables")
copy_files(src_cal, DEST_CALTABLES_LUSTRE, "Step 1b: Copying calibration tables to lustre/caltables")

# Step 2: Copy src_cal_beam/* to caltables_beam directories
copy_files(src_cal_beam, DEST_CALTABLES_BEAM, "Step 2a: Copying beam calibration tables to caltab/caltables_beam")
copy_files(src_cal_beam, DEST_CALTABLES_BEAM_LUSTRE, "Step 2b: Copying beam calibration tables to lustre/caltables_beam")

# Step 3: Remove old files from latest directories
remove_files(DEST_CALTABLES_LATEST, "Step 3a: Removing old files from caltab/caltables_latest")
remove_files(DEST_CALTABLES_BEAM_LATEST, "Step 3b: Removing old files from caltab/caltables_beam_latest")
remove_files(DEST_CALTABLES_LATEST_LUSTRE, "Step 3c: Removing old files from lustre/caltables_latest")
remove_files(DEST_CALTABLES_BEAM_LATEST_LUSTRE, "Step 3d: Removing old files from lustre/caltables_beam_latest")

# Step 4: Copy src_cal/* to latest directories
copy_files(src_cal, DEST_CALTABLES_LATEST, "Step 4a: Copying calibration tables to caltab/caltables_latest")
copy_files(src_cal, DEST_CALTABLES_LATEST_LUSTRE, "Step 4b: Copying calibration tables to lustre/caltables_latest")

# Step 5: Copy src_cal_beam/* to latest_beam directories
copy_files(src_cal_beam, DEST_CALTABLES_BEAM_LATEST, "Step 5a: Copying beam calibration tables to caltab/caltables_beam_latest")
copy_files(src_cal_beam, DEST_CALTABLES_BEAM_LATEST_LUSTRE, "Step 5b: Copying beam calibration tables to lustre/caltables_beam_latest")

# Step 6: Run pdsh commands for remote distribution
if not args.skip_pdsh:
    run_pdsh_command(
        "pdsh -w lwacalim[05-09] 'rm -r /fast/solarpipe/caltables/*.bcal'",
        "Step 6a: Removing old bcal files from remote nodes"
    )
    run_pdsh_command(
        "pdsh -w lwacalim[05-09] 'cp -r /lustre/solarpipe/realtime_pipeline/caltables_latest/* /fast/solarpipe/caltables/'",
        "Step 6b: Copying latest calibration tables to remote nodes"
    )
else:
    print("\nSkipping pdsh commands (--skip-pdsh flag set)")

print("\n" + "="*60)
print("Distribution complete!")
print("="*60)

