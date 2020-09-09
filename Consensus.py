#!/usr/bin/env python3
import sys
import collections
# This script reads a sdf file and makes the consensus between vina and Rfscore

sdffile = sys.argv[1]
outfile = sys.argv[2]
with open(outfile, 'w'):
    pass

rfscores = []
vinascores = []
plecscores = []
with open(sdffile, 'r') as f:
    for line in f:
        if '>  <rfscore_v3>' in line:
            rfscores.append(float(next(f).strip()))
        if 'VINA RESULT:' in line:
            vinascore = float(line.split()[2])
            if vinascore < 0:
                vinascores.append(vinascore)
        if '>  <PLEClinear_p5_l1_s65536>' in line:
            plecscores.append(float(next(f).strip()))


maxrf = max(rfscores)
minrf = min(rfscores)
rfrange = maxrf - minrf


maxvina = max(vinascores)
minvina = min(vinascores)
vinarange = maxvina - minvina

maxplec = max(plecscores)
minplec = min(plecscores)
plecrange = maxplec - minplec

def getnext(file):
    return next(file).strip()

molecules = set()
consensus_scores = {}
with open(sdffile, 'r') as f:
    outstring = []
    getmolname = True
    for line in f:
        if getmolname is True:
            getmolname = False
            molname = line.strip()
        if '>  <MODEL>' in line:
            model = getnext(f)
            outstring.append(line)
            outstring.append(str(model)+'\n')
            continue
        if '>  <rfscore_v3>' in line:
            rfscore = float(getnext(f))
            relative_rfscore = (rfscore - minrf) / rfrange
            outstring.append(line)
            outstring.append(str(rfscore)+'\n')
            outstring.append('\n')
            outstring.append('>  <relative_rfscore_v3>\n')
            outstring.append(str(relative_rfscore)+'\n')
            continue

        if '>  <PLEClinear_p5_l1_s65536>' in line:
            plecscore = float(getnext(f))
            relative_plecscore = (plecscore - minplec) / plecrange
            outstring.append(line)
            outstring.append(str(plecscore)+'\n')
            outstring.append('\n')
            outstring.append('>  <relative_plecscore>\n')
            outstring.append(str(relative_plecscore)+'\n')
            continue

            print('rf: '+str(relative_rfscore))
        if 'VINA RESULT:' in line:
            vinascore = float(line.split()[2])
            relative_vinascore = abs(vinascore - maxvina) / vinarange
            outstring.append(line)
            continue

        if '$$$$' not in line:
            outstring.append(line)
            continue
        if '$$$$' in line:
            consensus_vina_rf = (relative_rfscore + relative_vinascore) /2
            consensus_vina_rf_plec = (relative_rfscore + relative_vinascore + relative_plecscore) /3
            outstring[0] = molname.split('/')[0]+'_'+model+'\n'
            outstring.append('>  <vina_result>\n')
            outstring.append(str(vinascore)+'\n')
            outstring.append('\n')
            outstring.append('>  <relative_vina_result>\n')
            outstring.append(str(relative_vinascore)+'\n')
            outstring.append('\n')
            outstring.append('>  <consensus_vina_rf>\n')
            outstring.append(str(consensus_vina_rf)+'\n')
            outstring.append('\n')
            outstring.append('>  <consensus_vina_rf_plec>\n')
            outstring.append(str(consensus_vina_rf_plec)+'\n')
            outstring.append('\n')
            outstring.append('>  <ZINC_ID>\n')
            outstring.append(molname.split('_')[0]+'\n')
            outstring.append('\n')
            if molname.startswith('decoy'):
                outstring.append('>  <type>\n')
                outstring.append('decoy\n')
                outstring.append('\n')
            elif molname.startswith('active'):
                outstring.append('>  <type>\n')
                outstring.append('active\n')
                outstring.append('\n')
            outstring.append(line)
            if vinascore < 0:
                with open(outfile, 'a') as outf:
                    outf.write(''.join(outstring))
                    consensus_scores[outstring[0]] = consensus_vina_rf_plec
                    molecules.add(molname.split('_')[0])
            outstring = []
            getmolname = True
            continue

if len(sys.argv) > 3:
    with open('top'+str(sys.argv[3])+'percent.sdf', 'w'):
        pass
    nmol = int(sys.argv[3]) * len(molecules) / 100
    n=1
    print(len(molecules))
    presented_mols = set()
    topmols = set()
    oredered_poses = collections.OrderedDict(sorted(consensus_scores.items(), key=lambda x: x[1], reverse=True))
    for pose, score in oredered_poses.items():
        if n > nmol:
            break
        if pose.strip().split('_')[0] in presented_mols:
            continue
        else:
            print(str(n)+': '+pose.strip())
            topmols.add(pose.strip())
            presented_mols.add(pose.strip().split('_')[0])
            n += 1

    with open(outfile, 'r') as f:
        outstring = []
        getmolname = True
        for line in f:
            if getmolname is True:
                getmolname = False
                molname = line.strip()
            if '$$$$' not in line:
                outstring.append(line)
                continue
            if '$$$$' in line:
                outstring.append(line)
                if molname in topmols:
                    with open('top'+str(sys.argv[3])+'percent.sdf', 'a') as outf:
                        outf.write(''.join(outstring))
                outstring = []
                getmolname = True
                continue
