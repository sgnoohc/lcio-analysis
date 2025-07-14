import ROOT

ROOT.gROOT.SetBatch(True)

# Open the ROOT file
f = ROOT.TFile.Open("trackingNtuple.root")

# Variables to process
variables = ["pt", "eta", "phi"]

# Loop over variables
graphs = {}
c = ROOT.TCanvas("", "")
for var in variables:
    h_num = f.Get(f"num_{var}")
    h_den = f.Get(f"den_{var}")
    if not h_num or not h_den:
        print(f"Missing histogram for {var}, skipping")
        continue

    g_eff = ROOT.TGraphAsymmErrors(h_num, h_den)
    g_eff.SetTitle(f"Efficiency vs {var};{var};Efficiency")
    g_eff.SetMarkerStyle(20)
    g_eff.SetMarkerColor(ROOT.kBlue)
    g_eff.SetLineColor(ROOT.kBlue)
    # graphs[var] = g_eff
    xmin = h_den.GetXaxis().GetXmin()
    xmax = h_den.GetXaxis().GetXmax()
    g_eff.GetXaxis().SetLimits(xmin, xmax)
    g_eff.GetYaxis().SetLimits(0, 1.2)

    # Draw and save each plot
    if var == "pt":
        c.SetLogx(True)
    else:
        c.SetLogx(False)
    g_eff.Draw("AP")
    c.SaveAs(f"eff_{var}.pdf")

