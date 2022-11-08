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

phylum = input(colored('Is your organism an [animal], [plant], [fungus], [protist], or [metazoan]? Input one of the options within the brackets: ', 'cyan'))
phylum = phylum.lower()
organism_name = input(colored("Input the 'Genus species' for your reference organism: ", 'magenta'))
org_split = organism_name.lower().split()
org_dir = "_".join(org_split)

# This is where the user-friendliness comes in: use Ensembl's highly regular patterning (for eukaryotes) to fetch the needed file for supported organisms.

# Vertebrates/Model Animals:
if 'ani' in phylum:
   try:
       with ftplib.FTP('ftp.ensembl.org') as ftp:
           ftp.login('anonymous')
           ftp.cwd('/pub/current_fasta/{}/cdna/'.format(org_dir))
           for filename in ftp.nlst(pattern):
               fhandle = open(filename, 'wb')
               ftp.retrbinary('RETR ' + filename, fhandle.write)
               fhandle.close()
   except ftplib.error_perm:
       text = colored("It appears that your organism is not supported. Try running again with a closely related species or check spelling.", 'red')
       print(text)
       sys.exit(1)

# Fungi:
if 'fu' in phylum:
   try:
       with ftplib.FTP('ftp.ensemblgenomes.org') as ftp:
           ftp.login('anonymous')
           ftp.cwd('/pub/current/fungi/fasta/{}/cdna/'.format(org_dir))
           for filename in ftp.nlst(pattern):
               fhandle = open(filename, 'wb')
               ftp.retrbinary('RETR ' + filename, fhandle.write)
               fhandle.close()
   except ftplib.error_perm:
       text = colored("It appears that your organism is not supported. Try running again with a closely related species or check spelling.", 'red')
       print(text)
       sys.exit(1)

# Plants:
if 'pl' in phylum:
   try:
       with ftplib.FTP('ftp.ensemblgenomes.org') as ftp:
           ftp.login('anonymous')
           ftp.cwd('/pub/current/plants/fasta/{}/cdna/'.format(org_dir))
           for filename in ftp.nlst(pattern):
               fhandle = open(filename, 'wb')
               ftp.retrbinary('RETR ' + filename, fhandle.write)
               fhandle.close()
   except ftplib.error_perm:
       text = colored("It appears that your organism is not supported. Try running again with a closely related species or check spelling.", 'red')
       print(text)
       sys.exit(1)

# Protists:
if 'pr' in phylum:
    try:
        with ftplib.FTP('ftp.ensemblgenomes.org') as ftp:
            ftp.login('anonymous')
            ftp.cwd('/pub/current/protists/fasta/{}/cdna/'.format(org_dir))
            for filename in ftp.nlst(pattern):
                fhandle = open(filename, 'wb')
                ftp.retrbinary('RETR ' + filename, fhandle.write)
                fhandle.close()
    except ftplib.error_perm:
        text = colored("It appears that your organism is not supported. Try running again with a closely related species or check spelling.", 'red')
        print(text)
        sys.exit(1)

# Metazoans:
if 'met' in phylum:
    try:
        with ftplib.FTP('ftp.ensemblgenomes.org') as ftp:
            ftp.login('anonymous')
            ftp.cwd('/pub/current/metazoa/fasta/{}/cdna/'.format(org_dir))
            for filename in ftp.nlst(pattern):
                fhandle = open(filename, 'wb')
                ftp.retrbinary('RETR ' + filename, fhandle.write)
                fhandle.close()
    except ftplib.error_perm:
        text = colored("It appears that your organism is not supported. Try running again with a closely related species or check spelling.", 'red')
        print(text)
        sys.exit(1)

# Error Handling:
valid_responses = ('ani', 'fu', 'pl')
if not any(s in phylum for s in valid_responses):
   text = colored("Looks like you've entered an improper phylum. Try again!", "red")
   print(text)
   sys.exit(1)

# Looks like everything works up to here now.

# Define some potentially useful variables that I can plug into subprocess.
kallisto = './kallisto'
index = 'index'
quant = 'quant'

# Fetch the index file (name independent, it just needs to end in the pattern defined above)
# All Ensembl cDNA files should have that precise formatting.

pattern = '*.cdna.all.fa.gz'

ref_cdna = glob.glob(pattern)
ref_cdna = ref_cdna[0]

exts = ('*.fastq.gz', '*.fq.gz', "*.fastq", "*.fq")

threads = int(os.cpu_count())

# This checks whether you already have an index and, if so, whether you want to build a new one.
index_checker = os.path.isfile('kallisto_index')

if index_checker:
    indexed = input(colored("Would you like to build a fresh transcriptome index? This is optional [Y/N] ", 'green'))
    if indexed == "Y":
        subprocess.call([kallisto, index, '-i', 'kallisto_index', ref_cdna])
    else:
        pass
else:
    subprocess.call([kallisto, index, '-i', 'kallisto_index', ref_cdna])

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
    fastqs = glob.glob("./*_1.f*q*")
    try:
        for forward in fastqs:
            reverse = forward.replace("_1.f", "_2.f")
            matcher = re.search("_1.f.*q*", forward).group(0)
            dir = forward.replace(matcher, "")
            subprocess.run([kallisto, quant,
            '-i', 'kallisto_index',
            '-o', dir + '_quant',
            '--bias',
            '-b', '200',
            '-t', threads,
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
    fastqs = glob.glob("./*.f*q*")
    frag_len = 150
    standard_dev = 10
    frags = input("If known, input estimated average fragment length (from FastQC). If not known, hit enter: ")
    if not frags:
        frags = frag_len
    sd = input("If known, input the standard deviation of the fragment length (from FastQC). If not known, hit enter: ")
    if not sd:
        sd = standard_dev
    try:
        for read in fastqs:
            matcher = re.search(".f.*q*", read).group(0)
            dir = read.replace(matcher, "")
            subprocess.run([kallisto, quant,
            '-i', 'kallisto_index',
            '-o', dir + '_quant',
            '--bias',
            '-b', '200',
            '-t', threads,
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
