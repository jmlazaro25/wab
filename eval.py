import os
import sys
import numpy as np
import pickle as pkl
import xgboost as xgb
import mods.ROOTmanager as manager
from treeMaker import branches_info

pkl_file   = os.getcwd()+'/gabe.pkl'
model = pkl.load(open(pkl_file,'rb'))

def main():

    # Inputs and their trees and stuff
    pdict = manager.parse()
    inlist = pdict['inlist']
    outlist = pdict['outlist']
    group_labels = pdict['groupls']
    maxEvent = pdict['maxEvent']

    branches_info['discValue_EcalVeto'] = {'rtype': float, 'default': 0.5}

    # Construct tree processes
    procs = []
    for gl, group in zip(group_labels, inlist):
        procs.append( manager.TreeProcess(event_process, group, ID=gl, tree_name='EcalVeto',
            pfreq=100) )

    # Process jobs
    for proc in procs:

        print('\nRunning %s'%(proc.ID))
        
        # Move into appropriate scratch dir
        os.chdir(proc.tmp_dir)

        # Make an output file and new tree (copied from input + discValue)
        proc.tfMaker = manager.TreeMaker(group_labels[procs.index(proc)]+'.root',\
                                         "EcalVeto",\
                                         branches_info,\
                                         outlist[procs.index(proc)]
                                         )

        # RUN
        proc.extrafs = [ proc.tfMaker.wq ] # Gets executed at the end of run()
        proc.run(maxEvents=maxEvent)

    # Remove scratch directory if there is one
    manager.rmScratch()

    print('\nDone!\n')


def event_process(self):

    # Feature list from input tree
    # Exp: feats = [ feat_value for feat_value in self.tree~ ]
    # Put all segmentation variables in for now (Take out the ones we won't need once
    # we make sure that all the python bdt stuff works)
    feats = [
            # Base variables
            self.tree.nReadoutHits              ,
            self.tree.summedDet                 ,
            self.tree.summedTightIso            ,
            self.tree.maxCellDep                ,
            self.tree.showerRMS                 ,
            self.tree.xStd                      ,
            self.tree.yStd                      ,
            self.tree.avgLayerHit               ,
            self.tree.stdLayerHit               ,
            self.tree.deepestLayerHit           ,
            self.tree.ecalBackEnergy            ,
            # Electron RoC variables
            self.tree.eContEnergy_x1_s1         ,
            self.tree.eContEnergy_x2_s1         ,
            self.tree.eContEnergy_x3_s1         ,
            self.tree.eContEnergy_x4_s1         ,
            self.tree.eContEnergy_x5_s1         ,
            # Photon RoC variables
            self.tree.gContEnergy_x1_s1         ,
            self.tree.gContEnergy_x2_s1         ,
            self.tree.gContEnergy_x3_s1         ,
            self.tree.gContEnergy_x4_s1         ,
            self.tree.gContEnergy_x5_s1         ,
            # Outside RoC variables
            self.tree.oContEnergy_x1_s1         ,
            self.tree.oContEnergy_x2_s1         ,
            self.tree.oContEnergy_x3_s1         ,
            self.tree.oContEnergy_x4_s1         ,
            self.tree.oContEnergy_x5_s1         ,
            self.tree.oContNHits_x1_s1          ,
            self.tree.oContNHits_x2_s1          ,
            self.tree.oContNHits_x3_s1          ,
            self.tree.oContNHits_x4_s1          ,
            self.tree.oContNHits_x5_s1          ,
            self.tree.oContXStd_x1_s1           ,
            self.tree.oContXStd_x2_s1           ,
            self.tree.oContXStd_x3_s1           ,
            self.tree.oContXStd_x4_s1           ,
            self.tree.oContXStd_x5_s1           ,
            self.tree.oContYStd_x1_s1           ,
            self.tree.oContYStd_x2_s1           ,
            self.tree.oContYStd_x3_s1           ,
            self.tree.oContYStd_x4_s1           ,
            self.tree.oContYStd_x5_s1
            ]

    # Copy input tree feats to new tree
    for feat_name, feat_value in zip(self.tfMaker.branches_info, feats):
        self.tfMaker.branches[feat_name][0] = feat_value

    # Add prediction to new tree
    evtarray = np.array([feats])
    pred = float(model.predict(xgb.DMatrix(evtarray))[0])
    self.tfMaker.branches['discValue_EcalVeto'][0] = pred

    # Fill new tree with current event values
    self.tfMaker.tree.Fill()

if __name__ == "__main__":
    main()
