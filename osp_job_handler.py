import argparse
import socket
import time

def main():
    parser = argparse.ArgumentParser(description="OSP Task Handler")
    parser.add_argument("task_id", type=int, help="Task ID")
    parser.add_argument("task_count", type=int, help="Total number of tasks")
    args = parser.parse_args()

    hostname = socket.gethostname()

    print(f"Task ID: {args.task_id}, Task Count: {args.task_count}, Hostname: {hostname}")
    
    # Sleep for 3 seconds
    time.sleep(3)

    # Add your task-specific code here
    # For example:
    # process_data(args.task_id, args.task_count)

if __name__ == "__main__":
    main()