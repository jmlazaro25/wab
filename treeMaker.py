import os
import math
import ROOT as r
import numpy as np
from mods import ROOTmanager as manager
from mods import physTools, mipTracking

cellMap = np.loadtxt( 'mods/cellmodule.txt' )

r.gSystem.Load(
    '/home/jmlazaro/research/ldmx-sw-v2.3.0/install/lib/libEvent.so'
    )

# TreeModel to build here
branches_info = {
    # Base variables
    'nReadoutHits': {
        'rtype': int,
        'default': 0
        },
    'summedDet': {
        'rtype': float,
        'default': 0.
        },
    'summedTightIso': {
        'rtype': float,
        'default': 0.
        },
    'maxCellDep': {
        'rtype': float,
        'default': 0.
        },
    'showerRMS': {
        'rtype': float,
        'default': 0.
        },
    'xStd': {
        'rtype': float,
        'default': 0.
        },
    'yStd': {
        'rtype': float,
        'default': 0.
        },
    'avgLayerHit': {
        'rtype': float,
        'default': 0.
        },
    'stdLayerHit': {
        'rtype': float,
        'default': 0.
        },
    'deepestLayerHit': {
        'rtype': int,
        'default': 0
        },
    'ecalBackEnergy': {
        'rtype': float,
        'default': 0.
        },
    # HcalVeto
    'maxPE': {
        'rtype': int,
        'default': 0.5
        },
    # Extra Momentum and Thetale Info for analysis only (not trained on)
    'recoilPT': {
        'rtype': float,
        'default': 0.
        },
    'eTheta': {
        'rtype': float,
        'default': 0.
        },
    'gTheta': {
        'rtype': float,
        'default': 0.
        },
    'egTheta': {
        'rtype': float,
        'default': 0.
        }
    }

for i in range( physTools.nRegions ):
    branches_info[ 'electronContainmentEnergy_x{}'.format( i + 1 ) ] = {
        'rtype': float,
        'default': 0.
        }
    branches_info[ 'photonContainmentEnergy_x{}'.format( i + 1 ) ] = {
        'rtype': float,
        'default': 0.
        }
    branches_info[ 'outsideContainmentEnergy_x{}'.format( i + 1 ) ] = {
        'rtype': float,
        'default': 0.
        }
    branches_info[ 'outsideContainmentNHits_x{}'.format( i + 1 ) ] = {
        'rtype': int,
        'default': 0
        }
    branches_info[ 'outsideContainmentXStd_x{}'.format( i + 1 ) ] = {
        'rtype': float,
        'default': 0.
        }
    branches_info[ 'outsideContainmentYStd_x{}'.format( i + 1 ) ] = {
        'rtype': float,
        'default': 0.
        }


def main():

    # Inputs and their trees and stuff
    pdict = manager.parse()
    batch_mode = pdict[ 'batch' ]
    separate = pdict[ 'separate' ]
    inlist = pdict[ 'inlist' ]
    outlist = pdict[ 'outlist' ]
    group_labels = pdict[ 'groupls' ]
    maxEvents = pdict[ 'maxEvents' ]
    # Should maybe put in parsing eventually and make event_process *arg

    # Construct tree processes
    procs = []
    for gl, group in zip( group_labels, inlist ):
        procs.append(
            manager.TreeProcess(
                event_process, group, ID=gl, batch=batch_mode, pfreq=100
                )
            )

    # Process jobs
    for proc in procs:

        # Move into appropriate scratch dir
        os.chdir( proc.tmp_dir )

        # Branches needed
        proc.ecalVeto = proc.addBranch( 'EcalVetoResult', 'EcalVeto_v12' )
        proc.hcalVeto = proc.addBranch( 'HcalVetoResult', 'HcalVeto_v12' )
        proc.targetSPHits = proc.addBranch(
            'SimTrackerHit', 'TargetScoringPlaneHits_v12'
            )
        proc.ecalSPHits = proc.addBranch(
            'SimTrackerHit', 'EcalScoringPlaneHits_v12'
            )
        proc.ecalRecHits = proc.addBranch( 'EcalHit', 'EcalRecHits_v12' )

        # Tree/Files(s) to make
        print( '\nRunning %s' % ( proc.ID ) )

        proc.separate = separate

        proc.tfMakers = { 'unsorted': None}
        if proc.separate:
            proc.tfMakers = {
                'egin': None,
                'ein': None,
                'gin': None,
                'none': None
                }

        for tfMaker in proc.tfMakers:
            proc.tfMakers[tfMaker] = manager.TreeMaker(group_labels[procs.index(proc)]+\
                                        '_{}.root'.format(tfMaker),\
                                        "EcalVeto",\
                                        branches_info,\
                                        outlist[procs.index(proc)]
                                        )

        # Gets executed at the end of run()
        proc.extrafs = [
            proc.tfMakers[ tfMaker ].wq for tfMaker in proc.tfMakers
            ]

        # RUN
        proc.run( maxEvents=maxEvents )

    # Remove scratch directory if there is one
    if not batch_mode: # Don't want to break other batch jobs when one finishes
        manager.rmScratch()

    print( '\nDone!\n' )


# Process an event
def event_process( self ):

    # Initialize BDT input variables w/ defaults
    feats = next( iter( self.tfMakers.values() ) ).resetFeats()

    # Assign pre-computed variables
    feats[ 'nReadoutHits' ] = self.ecalVeto.getNReadoutHits()
    feats[ 'summedDet' ] = self.ecalVeto.getSummedDet()
    feats[ 'summedTightIso' ] = self.ecalVeto.getSummedTightIso()
    feats[ 'maxCellDep' ] = self.ecalVeto.getMaxCellDep()
    feats[ 'showerRMS' ] = self.ecalVeto.getShowerRMS()
    feats[ 'xStd' ] = self.ecalVeto.getXStd()
    feats[ 'yStd' ] = self.ecalVeto.getYStd()
    feats[ 'avgLayerHit' ] = self.ecalVeto.getAvgLayerHit()
    feats[ 'stdLayerHit' ] = self.ecalVeto.getStdLayerHit()
    feats[ 'deepestLayerHit' ] = self.ecalVeto.getDeepestLayerHit()
    feats[ 'ecalBackEnergy' ] = self.ecalVeto.getEcalBackEnergy()
    feats[ 'maxPE' ] = self.hcalVeto.getMaxPEHit().getPE()

    for i in range( physTools.nRegions ):
        feats[ 'electronContainmentEnergy_x{}'.format( i + 1 )
              ] = self.ecalVeto.getElectronContainmentEnergy()[ i ]
        feats[ 'photonContainmentEnergy_x{}'.format( i + 1 )
              ] = self.ecalVeto.getPhotonContainmentEnergy()[ i ]
        feats[ 'outsideContainmentEnergy_x{}'.format( i + 1 )
              ] = self.ecalVeto.getOutsideContainmentEnergy()[ i ]
        feats[ 'outsideContainmentNHits_x{}'.format( i + 1 )
              ] = self.ecalVeto.getOutsideContainmentNHits()[ i ]
        feats[ 'outsideContainmentXStd_x{}'.format( i + 1 )
              ] = self.ecalVeto.getOutsideContainmentXStd()[ i ]
        feats[ 'outsideContainmentYStd_x{}'.format( i + 1 )
              ] = self.ecalVeto.getOutsideContainmentYStd()[ i ]

    ###################################
    # Determine event type
    ###################################

    # Get e position and momentum from EcalSP
    e_ecalHit = physTools.electronEcalSPHit( self.ecalSPHits )
    if e_ecalHit != None:
        e_ecalPos, e_ecalP = e_ecalHit.getPosition(), e_ecalHit.getMomentum()

        # Extra analysis info
        feats[ 'recoilPT' ] = physTools.mag( e_ecalP )

    # Photon Info from targetSP
    e_targetHit = physTools.electronTargetSPHit( self.targetSPHits )
    if e_targetHit != None:
        g_targPos, g_targP = physTools.gammaTargetInfo( e_targetHit )

        # Extra analysis info
        feats[ 'eTheta' ] = physTools.angle(
            e_targetHit.getMomentum(), units='degrees'
            )
        feats[ 'gTheta' ] = physTools.angle( g_targP, units='degrees' )
        feats[ 'egTheta' ] = physTools.angle(
            g_targP, units='degrees', vec2=e_targetHit.getMomentum()
            )

    else: # Should about never happen -> division by 0 in g_traj
        print( 'no e at targ!' )
        g_targPos = g_targP = np.zeros( 3 )

    # Get electron and photon trajectories AND
    # Fiducial categories (filtered into different output trees)
    if self.separate:
        e_traj = g_traj = None
        e_fid = g_fid = False

        if e_ecalHit != None:
            e_traj = physTools.layerIntercepts( e_ecalPos, e_ecalP )
            for cell in cellMap:
                if physTools.dist(
                    cell[ 1: ], e_traj[ 0 ]
                    ) <= physTools.cell_radius:
                    e_fid = True
                    break

        if e_targetHit != None:
            g_traj = physTools.layerIntercepts( g_targPos, g_targP )
            for cell in cellMap:
                if physTools.dist(
                    cell[ 1: ], g_traj[ 0 ]
                    ) <= physTools.cell_radius:
                    g_fid = True
                    break

    ###################################
    # Fill the tree (according to fiducial category) with values for this event
    ###################################
    if not self.separate:
        self.tfMakers[ 'unsorted' ].fillEvent( feats )
    else:
        if e_fid and g_fid: self.tfMakers[ 'egin' ].fillEvent( feats )
        elif e_fid and not g_fid: self.tfMakers[ 'ein' ].fillEvent( feats )
        elif not e_fid and g_fid: self.tfMakers[ 'gin' ].fillEvent( feats )
        else: self.tfMakers[ 'none' ].fillEvent( feats )


if __name__ == "__main__":
    main()
