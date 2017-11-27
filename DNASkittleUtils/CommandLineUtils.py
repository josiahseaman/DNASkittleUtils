
from __future__ import print_function

import os
import subprocess
import sys
import datetime

output_dir = os.getcwd()
# one log file per day across all jobs
log_file_name = "WHAT_I_DID_" + str(datetime.date.today()) + '.log'


def log_command(args):
    command = ' '.join(args) if isinstance(args, list) else args
    print(command)
    with open(os.path.join(output_dir, log_file_name), 'a') as log:
        log.write(command + '\n')
    return command


def call_output(args):
    command = log_command(args)
    return subprocess.check_output(command, shell=True)


def call(args):
    command = log_command(args)
    return subprocess.call(command, shell=True)
    # TODO: add error handling for bad return code


def remove_extensions(path):
    """Remove extension and path"""
    first_extension = os.path.splitext(path)[0]
    return os.path.splitext(first_extension)[0]


def just_the_name(path):
    """Remove extension and path"""
    return remove_extensions(os.path.basename(path))


def delete_file_contents(file_path, scratch_only=False):
    """When the presence of a file is being used as an indicator of what files have already been computed,
    we want to keep the file even after it has already been used in the next processing step.  This
    function deletes the contents of the file while leaving it in its place.
    scratch_only: bool whether file deletion 'file_path' must contain "scratch" in the name
    """
    if not scratch_only or 'scratch' in file_path:
        if os.path.exists(file_path):
            with open(file_path, 'w') as big_file:
                big_file.write('Contents deleted to save scratch space')
                print('File contents deleted:', file_path)
    elif scratch_only:
        print("ERROR: Not blanking file because it's not in a scratch folder", file_path, file=sys.stderr)


def make_output_dir_with_suffix(base_path, suffix):
    output_dir = base_path + suffix
    print("Creating Directory...", os.path.basename(output_dir))
    os.makedirs(output_dir, exist_ok=True)
    return output_dir