import os
import sys
import math
import numpy as np
import pickle as pkl
import xgboost as xgb
import mods.ROOTmanager as manager
from gabTreeMaker import branches_info

pkl_file = os.getcwd() + '/gabrielle_p.pkl'
model = pkl.load(open(pkl_file, 'rb'))


def main():

    # Inputs and their trees and stuff
    pdict = manager.parse()
    inlist = pdict['inlist']
    outlist = pdict['outlist']
    group_labels = pdict['groupls']
    maxEvents = pdict['maxEvents']

    # Add discValue to tree model
    branches_info['discValue_EcalVeto'] = {'rtype': float, 'default': 0.5}
    # Construct tree processes
    procs = []
    for gl, group in zip(group_labels, inlist):
        procs.append(
            manager.TreeProcess(event_process,
                                group,
                                ID=gl,
                                tree_name='EcalVeto',
                                pfreq=100))

    # Process jobs
    for proc in procs:

        print('\nRunning %s' % (proc.ID))

        # Move into appropriate scratch dir
        os.chdir(proc.tmp_dir)

        # Make an output file and new tree (copied from input + discValue)
        proc.tfMaker = manager.TreeMaker(group_labels[procs.index(proc)]+'.root',\
                                         "EcalVeto",\
                                         branches_info,\
                                         outlist[procs.index(proc)]
                                         )

        # RUN
        proc.extrafs = [proc.tfMaker.wq]  # Gets executed at the end of run()
        proc.run(maxEvents=maxEvents)

    # Remove scratch directory if there is one
    manager.rmScratch()

    print('\nDone!\n')


def event_process(self):

    # Feature list from input tree
    # Exp: feats = [ feat_value for feat_value in self.tree~ ]
    feats = [
        # Base variables
        self.tree.nReadoutHits,  # 0
        self.tree.summedDet,  # 1
        self.tree.summedTightIso,  # 2
        self.tree.maxCellDep,  # 3
        self.tree.showerRMS,  # 4
        self.tree.xStd,  # 5
        self.tree.yStd,  # 6
        self.tree.avgLayerHit,  # 7
        self.tree.stdLayerHit,  # 8
        self.tree.deepestLayerHit,  # 9
        self.tree.ecalBackEnergy,  # 10
        # Electron RoC variables
        self.tree.electronContainmentEnergy_x1,  # 11
        self.tree.electronContainmentEnergy_x2,  # 12
        self.tree.electronContainmentEnergy_x3,  # 13
        self.tree.electronContainmentEnergy_x4,  # 14
        self.tree.electronContainmentEnergy_x5,  # 15
        # Photon RoC variables
        self.tree.photonContainmentEnergy_x1,  # 16
        self.tree.photonContainmentEnergy_x2,  # 17
        self.tree.photonContainmentEnergy_x3,  # 18
        self.tree.photonContainmentEnergy_x4,  # 19
        self.tree.photonContainmentEnergy_x5,  # 20
        # Outside RoC variables
        self.tree.outsideContainmentEnergy_x1,  # 21
        self.tree.outsideContainmentEnergy_x2,  # 22
        self.tree.outsideContainmentEnergy_x3,  # 23
        self.tree.outsideContainmentEnergy_x4,  # 24
        self.tree.outsideContainmentEnergy_x5,  # 25
        self.tree.outsideContainmentNHits_x1,  # 26
        self.tree.outsideContainmentNHits_x2,  # 27
        self.tree.outsideContainmentNHits_x3,  # 28
        self.tree.outsideContainmentNHits_x4,  # 29
        self.tree.outsideContainmentNHits_x5,  # 30
        self.tree.outsideContainmentXStd_x1,  # 31
        self.tree.outsideContainmentXStd_x2,  # 32
        self.tree.outsideContainmentXStd_x3,  # 33
        self.tree.outsideContainmentXStd_x4,  # 34
        self.tree.outsideContainmentXStd_x5,  # 35
        self.tree.outsideContainmentYStd_x1,  # 36
        self.tree.outsideContainmentYStd_x2,  # 37
        self.tree.outsideContainmentYStd_x3,  # 38
        self.tree.outsideContainmentYStd_x4,  # 39
        self.tree.outsideContainmentYStd_x5  # 40
    ]

    # Copy input tree feats to new tree
    for feat_name, feat_value in zip(self.tfMaker.branches_info, feats):

        # Condition to catch NaN bugs if they exist
        if math.isnan(feat_value):
            print('The problematic feature is: {}'.format(feat_name))
            self.tfMaker.branches[feat_name][0] = self.tfMaker.branches_info[
                feat_name]['default']
        else:
            self.tfMaker.branches[feat_name][0] = feat_value

    # Add the quantities needed for BDT analysis to the tree
    self.tfMaker.branches['maxPE'][0] = self.tree.maxPE
    self.tfMaker.branches['recoilPT'][0] = self.tree.recoilPT
    self.tfMaker.branches['eTheta'][0] = self.tree.eTheta
    self.tfMaker.branches['gTheta'][0] = self.tree.gTheta
    self.tfMaker.branches['egTheta'][0] = self.tree.egTheta

    # Add prediction to new tree
    evtarray = np.array([feats])
    pred = float(model.predict(xgb.DMatrix(evtarray))[0])
    self.tfMaker.branches['discValue_EcalVeto'][0] = pred

    # Fill new tree with current event values
    self.tfMaker.tree.Fill()


if __name__ == "__main__":
    main()
