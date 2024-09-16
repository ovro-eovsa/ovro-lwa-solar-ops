import os
import socket
import time

def main():
    # Get the task ID and total number of tasks from SLURM environment variables
    task_id = int(os.environ.get('SLURM_PROCID', 0))
    task_count = int(os.environ.get('SLURM_NTASKS', 1))

    hostname = socket.gethostname()
<<<<<<<< HEAD:slurm_taskid_test.py
========

    print(f"Task ID: {task_id}, Task Count: {task_count}, Hostname: {hostname}")
>>>>>>>> 5cb2ba12ccb9378d5565d270a141596ceaa425b9:osp_task_handler.py
    
    print(time.ctime(),f"SLURM_PROCID: {task_id}, SLURM_NTASKS: {task_count}, Hostname: {hostname}")
    # should print correct task_id and task_count and hostname

    # Sleep for 3 seconds
    time.sleep(3)

<<<<<<<< HEAD:slurm_taskid_test.py
========
    # Add your task-specific code here
    # For example:
    # process_data(task_id, task_count)

>>>>>>>> 5cb2ba12ccb9378d5565d270a141596ceaa425b9:osp_task_handler.py
if __name__ == "__main__":
    main()