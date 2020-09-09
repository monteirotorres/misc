#!/usr/bin/env python

import multiprocessing as mp
import subprocess
import os
global receptor


def rescore(folder):
    try:
        pdbqtpath = os.path.join(folder, "out.pdbqt")
        sdfpath = os.path.join(folder, "out.sdf")
        command = ["oddt_cli", "-i", "pdbqt", "--score", "rfscore_v3_pdbbind2016", "--score", "pleclinear_pdbbind2016", "--receptor", receptor, "-o", "sdf", "-O", sdfpath, pdbqtpath]
        subprocess.run(command)

    except subprocess.CalledProcessError:
        print('Command ' + ' '.join(command) + ' returned non-zero status\n' )
        pass


receptor = input("Please state the name of the receptor: ")

folderlist = []
for folder in os.listdir():
    if os.path.isfile(os.path.join(folder, "out.pdbqt")):
        folderlist.append(os.path.basename(folder))

p = mp.Pool(16)

for result in p.imap_unordered(rescore, folderlist, chunksize=100):
    pass
