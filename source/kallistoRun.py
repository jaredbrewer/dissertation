#!/usr/bin/env python3

# macOS ONLY

# People need to have installed Anaconda for this to work. Presumably that can be done easily enough. Alt: execute another script first.

import re, os, ftplib, subprocess, glob, sys, shutil, platform

# This gives the script some self awareness. It finds itself and changes the working directory to that path (temporarily).
# This is important for executing the brew_installer.sh script.
script_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(script_path)

# Welcome to the era of Apple Silicon.

if "macOS" in platform.platform():
    # This will check to see if several important system programs are installed in sequential order.
    # If they are not, then it executes a script to install them, may require user password.
    if not shutil.which('xcode-select'):
        subprocess.run(['./brew_installer.sh'])
    if not shutil.which('brew'):
        subprocess.run(['./brew_installer.sh'])
    else:
        subprocess.run(["brew", "install", "hdf5"])
    if "arm64" in platform.platform():
        kallisto = "./kallisto_as"
    elif "x86_64" in platform.platform():
        kallisto = "./kallisto_intel"
elif "Linux" in platform.platform():
    kallisto = "./kallisto_linux"
    subprocess.run(["sudo", "apt-get", "install", "libhdf5-dev"])

# This will check to see if several important system programs are installed in hierarchical order
# If they are not, then it executes a script to install them, may require user password.

from termcolor import colored, cprint
import Bio

fastq_dir = input("Enter the directory of your FASTQ files (drag and drop is fine): ")
fastq_dir = fastq_dir.strip()
try:
    os.chdir(fastq_dir)
except FileNotFoundError:
    text = colored("Looks like that directory does not exist - restart the script and try dragging and dropping the folder directly into Terminal.", "red")
    print(text)
    sys.exit(1)
except PermissionError:
    print(text)
    sys.exit(1)

# Define some potentially useful variables that I can plug into subprocess.
# A version of Kallisto built with HDF5 is essential - I am a little bit lost on why they stopped bundling HDF5 with newer version with no bootstrapping replacement. A version of Kallisto with HDF5 for macOS is included, but this is not a sustainable distribution method.

index = 'index'
quant = 'quant'

# This checks whether you already have an index and, if so, whether you want to build a new one.
index_checker = os.path.isfile('kallisto_index')

if index_checker:
    indexed = input(colored("Would you like to build a fresh transcriptome index? This is optional [Y/N] ", 'green'))
    if indexed == "Y":
        ref_cdna = input(colored("Provide a path to the reference cDNA for your organism: ", "green")).rstrip()
        # This defies logic - why does the drag-and-drop add a space? I guess no harm in stripping either way, but so weird as far as a source of error.
        subprocess.run([kallisto, index, '-i', 'kallisto_index', ref_cdna])
    else:
        pass
else:
    ref_cdna = input(colored("Provide a path to the reference cDNA for your organism: ", "green")).rstrip()
    subprocess.run([kallisto, index, '-i', 'kallisto_index', ref_cdna])

while True:
    valid = ('pe', 'se')
    read_type = input(colored("Are your reads paired-end [PE] or single-end [SE]? Paired-end reads have two files for each condition - one with '_1' and one with '_2' [PE/SE] ", "yellow"))
    if read_type.lower() not in valid:
        text = colored("Looks like you had a typo in the last prompt!", "red")
        print(text)
        continue
    else:
        break

# This script runs for PE reads - it can estimate fragment length from this and doesn't need you to provide the information.

cnt = 0

if 'PE' in read_type.upper():
    fastqs = glob.glob("./*_1.f*q.gz")
    try:
        for forward in fastqs:
            reverse = forward.replace("_1.f", "_2.f")
            matcher = re.search("_1.f.*q.gz", forward).group(0)
            dir = forward.replace(matcher, "")
            subprocess.run([kallisto, quant,
            '-i', 'kallisto_index',
            '-o', dir + '_quant',
            '--bias',
            '-b', '200',
            '-t', '4',
            forward, reverse])
            cnt += 1
            print(cnt)
    except:
        text = colored("Looks like something went wrong.", "red")
        print(text)
        sys.exit(1)
else:
    pass

# This needs a little bit more user-engagement - average fragment length and the SD are required, but can be substituted by guesses if it is not known (it rarely is).

cnt = 0

if 'SE' in read_type.upper():
    fastqs = glob.glob("./*.f*q.gz")
    frag_len = 150
    standard_dev = 10
    frags = input("If known, input estimated average fragment length (from FastQC). If not known, hit enter: ")
    if not frags:
        frags = frag_len
    sd = input("If known, input the standard deviation of the fragment length (from FastQC). If not known, hit enter: ")
    if not standard_dev:
        sd = standard_dev
    try:
        for read in fastqs:
            matcher = re.search(".f.*q.gz", read).group(0)
            dir = read.replace(matcher, "")
            subprocess.run([kallisto, quant,
            '-i', 'kallisto_index',
            '-o', dir + '_quant',
            '--bias',
            '-b', '200',
            '-t', '4',
            '--single',
            '-l', frags,
            '-s', sd,
            read])
            cnt += 1
            print(cnt)
    except:
        text = colored("Looks like something went wrong.", "red")
        print(text)
        sys.exit(1)
else:
    pass

input(colored("All finished! You are ready to proceed with further analysis and visualization in RStudio.", 'green', 'on_white'))
