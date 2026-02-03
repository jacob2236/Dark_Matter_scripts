import ROOT

f = ROOT.TFile.Open("reco.root")
#f.ls()  # lists contents

tree = f.Get("events")

for event in tree:
    particles = event.ReconstructedParticles
    for p in particles:
        if p.PDG == 11:   # electron
            mom = p.momentum
            px, py, pz = mom.x, mom.y, mom.z
            print("Electron p =", (px**2 + py**2 + pz**2)**0.5)