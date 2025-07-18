import ROOT

ROOT.gROOT.SetBatch(True)

# Open ROOT file
f = ROOT.TFile.Open("eff.root")

# Define system/layer structure from your map
nmodules_map = {
    1: {  # VB
        0: [16, 5, 1],
        1: [16, 5, 1],
        2: [15, 5, 1],
        3: [15, 5, 1],
        4: [21, 5, 1],
        5: [21, 5, 1],
        6: [29, 5, 1],
        7: [29, 5, 1],
    },
    2: {  # VE
        0: [1, 16, 2],
        1: [1, 16, 2],
        2: [1, 16, 2],
        3: [1, 16, 2],
        4: [1, 16, 2],
        5: [1, 16, 2],
        6: [1, 16, 2],
        7: [1, 16, 2],
    },
}

c = ROOT.TCanvas("c", "", 800, 600)

for system in nmodules_map:
    for layer in nmodules_map[system]:
        hnum_name = f"h_num_pt_sys{system}_layer{layer}"
        hden_name = f"h_den_pt_sys{system}_layer{layer}"
        h_num = f.Get(hnum_name)
        h_den = f.Get(hden_name)

        if not h_num or not h_den:
            print(f"Missing histograms for sys={system}, layer={layer}")
            continue

        g_eff = ROOT.TGraphAsymmErrors(h_num, h_den)
        g_eff.SetTitle(f"Efficiency vs pT (sys={system}, layer={layer});pT [GeV];Efficiency")
        g_eff.SetMarkerStyle(20)
        g_eff.SetMarkerColor(ROOT.kBlue)
        g_eff.SetLineColor(ROOT.kBlue)

        # X and Y axis limits
        xmin = h_den.GetXaxis().GetXmin()
        xmax = h_den.GetXaxis().GetXmax()
        g_eff.GetXaxis().SetLimits(xmin, xmax)
        g_eff.GetHistogram().SetMaximum(1.1)
        g_eff.GetHistogram().SetMinimum(0)

        # Draw
        c.SetLogx(True)
        g_eff.Draw("AP")
        c.SaveAs(f"plots/eff_pt_sys{system}_layer{layer}.pdf")
        c.SaveAs(f"plots/eff_pt_sys{system}_layer{layer}.png")

for system in nmodules_map:
    for layer in nmodules_map[system]:
        hnum_name = f"h_num_eta_sys{system}_layer{layer}"
        hden_name = f"h_den_eta_sys{system}_layer{layer}"
        h_num = f.Get(hnum_name)
        h_den = f.Get(hden_name)

        if not h_num or not h_den:
            print(f"Missing histograms for sys={system}, layer={layer}")
            continue

        g_eff = ROOT.TGraphAsymmErrors(h_num, h_den)
        g_eff.SetTitle(f"Efficiency vs pT (sys={system}, layer={layer});pT [GeV];Efficiency")
        g_eff.SetMarkerStyle(20)
        g_eff.SetMarkerColor(ROOT.kBlue)
        g_eff.SetLineColor(ROOT.kBlue)

        # X and Y axis limits
        xmin = h_den.GetXaxis().GetXmin()
        xmax = h_den.GetXaxis().GetXmax()
        g_eff.GetXaxis().SetLimits(xmin, xmax)
        g_eff.GetHistogram().SetMaximum(1.1)
        g_eff.GetHistogram().SetMinimum(0)

        # Draw
        c.SetLogx(False)
        g_eff.Draw("AP")
        c.SaveAs(f"plots/eff_eta_sys{system}_layer{layer}.pdf")
        c.SaveAs(f"plots/eff_eta_sys{system}_layer{layer}.png")

