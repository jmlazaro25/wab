import ROOT as r
import numpy as np
from array import array
import mods.ROOTmanager as manager
#from loguru import logger

r.gSystem.Load( '/home/jmlazaro/research/ldmx-sw-v2.3.0/install/lib/libEvent.so' )

#@logger.catch
def main():

    # Inputs and their trees and stuff
    pdict = manager.parse()
    inlist = pdict[ 'inlist' ]
    outlist = pdict[ 'outlist' ]
    group_labels = pdict[ 'groupls' ]
    maxEvents = pdict[ 'maxEvents' ]
    x_title = 'Maximum photoelectrons in an HCal Hit'
    y_title = 'Gabrielle discriminator value'
    # present as condition in https://github.com/LDMX-Software/Hcal/blob/9f968153cfe5683c94682c49f09451cce3b8cd25/src/Hcal/HcalVetoProcessor.cxx#L60-L72
    print( 'input group: {}'.format( inlist[0] ) ) # remove

    tree = manager.load(inlist[0], treeName='EcalVeto')

    # Construct tree processes
    procs = []
    for gl, group in zip( group_labels, inlist ):
        procs.append(
            manager.TreeProcess(
                event_process, group, ID=gl, tree = tree
                )
            )

    # Histograms, bramches, and stuff
    for proc in procs:

        # Misc
        print( 'Running %s' % ( proc.ID ) )
        proc.events_used = 0 # Left as a reminder that we might want to exclude some events
        proc.hists = {}
        #print(proc.tree.Print()) # remove

        # Histos
        bdtVhcal = manager.Histogram(r.TH2D(proc.ID,\
                ' ;' + x_title + ';' + y_title,
                20,0,100 , 20,0.95,1 ))
                #100,0,100 , 20,0.95,1 ))
                #100,0,100 , 100,0,100 ))
        proc.hists[ 'bdtVhcal' ] = bdtVhcal

        # RUN
        proc.run( maxEvents=maxEvents )

    # Gold
    plot( procs, outlist[ 0 ] + '/' + group_labels [ 0 ], x_title, y_title )

    manager.rmScratch() 

    print( '\nDone!\n' )

#@logger.catch
def event_process( self ):

    # Collect data
    self.hists[ 'bdtVhcal' ].hist.Fill(
        self.tree.maxPE, self.tree.discValue_gabrielle
        #self.tree.maxPE, self.tree.maxPE
        )


def plot( processes, output, xTitle='x', yTitle='y' ):

    c = r.TCanvas( "c", "", 900, 900 )
    c.SetTitle( ' ;' + xTitle + ';' + yTitle )
    c.SetLeftMargin( 0.15 )
    c.SetRightMargin( 0.15 )

    # Stylize and draw
    for proc in processes:
        for key in proc.hists:
            r.gStyle.SetOptStat( 0 )
            c.SetLogz()
            proc.hists[ key ].hist.SetMinimum( 1 )
            proc.hists[ key ].hist.Draw( "colz" )
            #print('\nUsed %s events\n' % (proc.events_used) )

    # Overall style
    #r.gStyle.SetOptTitle(0)

    # Save as pdf and png
    c.SaveAs( output + '.pdf' )


if __name__ == "__main__":
    main()
