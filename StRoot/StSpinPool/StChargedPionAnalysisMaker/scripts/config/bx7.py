name    = __name__.split('.')[-1]
VERSION = '$Id: bx7.py,v 1.1 2009/04/07 16:34:37 kocolosk Exp $'[5:-2]

import ROOT

class_ = ROOT.TH1D

binning = {
    'nbinsx': 128,
    'xbins': (-0.5, 127.5)
}

props = {
    'SetXTitle': (name,)
}

branches = ('mVertices*', 'mBbcTimeBin', 'mRunId', 'mBx7')

def accept_event(event):
    vertex_cut = event.nVertices() > 0
    real_data_cut = isinstance(event, ROOT.StChargedPionEvent)
    bin = event.bbcTimeBin()/32 
    bbc_cut = bin in (7,8,9) or (event.runId() > 7000000 and bin==6)
    return vertex_cut and real_data_cut and bbc_cut

def analyze(event, **kw):
    yield (event.bx7(),)

