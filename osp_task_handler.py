import os
import socket
import time

def main():
    # Get the task ID and total number of tasks from SLURM environment variables
    task_id = int(os.environ.get('SLURM_PROCID', 0))
    task_count = int(os.environ.get('SLURM_NTASKS', 1))
    hostname = socket.gethostname()
    print(f"SLURM_PROCID: {task_id}, SLURM_NTASKS: {task_count}, Hostname: {hostname}")
    
    # Sleep for 3 seconds
    time.sleep(3)

    # Add your task-specific code here
    # For example:
    # process_data(task_id, task_count)

if __name__ == "__main__":
    main()