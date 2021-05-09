import ROOT as r
import numpy as np
from array import array
import mods.ROOTmanager as manager
r.gSystem.Load('/home/jmlazaro/research/ldmx-sw/install/lib/libFramework.so')

def main():

    # Inputs and their trees and stuff
    pdict = manager.parse()
    inlist = pdict['inlist']
    outlist = pdict['outlist']
    group_labels = pdict['groupls']
    maxEvent = pdict['maxEvent']
    x_title = 'Maximum photoelectrons in an HCal Hit'
    y_title = 'BDT discriminator value'
    # present as condition in https://github.com/LDMX-Software/Hcal/blob/9f968153cfe5683c94682c49f09451cce3b8cd25/src/Hcal/HcalVetoProcessor.cxx#L60-L72

    # Construct tree processes
    procs = []
    for gl, group in zip(group_labels, inlist):
        procs.append( manager.TreeProcess(event_process, group,
                                          ID=gl, tree_name='EcalVeto') )

    # Histograms, bramches, and stuff
    for proc in procs:

        # Misc
        print('Running %s'%(proc.ID))
        proc.events_used = 0 # Left as a reminder that we might want to exclude some events
        proc.hists = {}

        # Histos
        bdtVhcal = manager.Histogram(r.TH2D(proc.ID,\
                ' ;' + x_title + ';' + y_title,
                20,0,100 , 20,0.95,1 ))
        proc.hists['bdtVhcal'] = bdtVhcal

        # RUN
        proc.run(maxEvents=maxEvent)

    # Gold
    plot(procs, outlist[0], x_title, y_title)

    print('\nDone!\n')

def event_process(self):

    # Collect data
    self.hists['bdtVhcal'].hist.Fill( self.tree.maxPE, self.tree.discValue_EcalVeto  )

def plot(processes, output, xTitle='x',yTitle='y'):

    c=r.TCanvas("c","",900,900)
    c.SetTitle(' ;' + xTitle + ';' + yTitle)
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.15)

    # Stylize and draw
    for proc in processes:
        for key in proc.hists:
            r.gStyle.SetOptStat(0)
            c.SetLogz()
            proc.hists[key].hist.SetMinimum(1)
            proc.hists[key].hist.Draw("colz")
            #print('\nUsed %s events\n' % (proc.events_used) )

    # Overall style
    #r.gStyle.SetOptTitle(0)

    # Save as pdf and png
    c.SaveAs(output + 'bdtVhcal.pdf')

if __name__ == "__main__":
    main()
