import ROOT

#---------------- *Open the ROOT files and get the histogram from each file* ----------

fIn1 = ROOT.TFile.Open("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_tttxtx_wminus_fullyLeptonic_signal/Events/run_01_decayed_1/unweighted_events.root")

h1 = fIn1.Get("wboson1_mass")
h2 = fIn1.Get("wboson2_mass")
h3 = fIn1.Get("wboson3_mass")
h4 = fIn1.Get("wboson4_mass")
h5 = fIn1.Get("wboson5_mass")

h1.SetLineWidth(2)
h2.SetLineWidth(2)
h3.SetLineWidth(2)
h4.SetLineWidth(2)
h5.SetLineWidth(2)

#--------------------------- *Scale the histograms by the cross sections* ------------
"""
xsec = [6.869e-00]
h1.Scale(xsec[0] / h1.Integral())
h2.Scale(xsec[1] / h2.Integral())
h3.Scale(xsec[2] / h3.Integral())
"""
#--------------------------- *Create a new canvas and draw the histograms* ------------

canvas = ROOT.TCanvas("canvas", "canvas", 600, 500)
h1.SetLineColor(ROOT.kRed)
h2.SetLineColor(ROOT.kBlue)
h3.SetLineColor(ROOT.kGreen)
h4.SetLineColor(ROOT.kBlack)
h5.SetLineColor(ROOT.kOrange)
h1.SetStats(0)

h1.GetXaxis().SetRangeUser(0, 100)  # Set x-axis limits for hist1
h1.GetYaxis().SetRangeUser(0,2600)  # Set y-axis limits for hist1

h1.Draw("E")
h2.Draw("E same")
h3.Draw("E same")
h4.Draw("E same")
h5.Draw("E same")

#-------------------------------- *Add the title and axis labels* ---------------------

h1.SetTitle("")
h1.GetXaxis().SetTitle("W_{mass} [GeV]")
h1.GetYaxis().SetTitle("Events")

tex1 = ROOT.TLatex(0.153,0.92,"t#bar{t}t#bar{t}w- signal")
tex1.SetNDC()
tex1.SetTextAngle(0)
tex1.SetTextFont(42)
tex1.SetTextSize(0.04)
tex1.SetTextAlign(11)
tex1.Draw()

tex2 = ROOT.TLatex(0.730,0.92,"\sqrt{S} = {13TeV}")
tex2.SetNDC()
tex2.SetTextAngle(0)
tex2.SetTextFont(42)
tex2.SetTextAlign(11)
tex2.SetTextSize(0.04)
tex2.Draw()


#--------------------------------------- *Add a legend* -------------------------------

legend = ROOT.TLegend(0.6, 0.7, 0.9, 0.9)
legend.AddEntry(h1, "W1 mass", "l")
legend.AddEntry(h2, "W2 mass", "l")
legend.AddEntry(h3, "W3 mass", "l")
legend.AddEntry(h4, "W4 mass", "l")
legend.AddEntry(h5, "W5 mass", "l")
legend.Draw()

#-------------------------------------- *Show the canvas* -----------------------------
#canvas.SetLogy()
canvas.Draw()
input("press enter to exit....")

#-------------------------------------------- *End* -----------------------------------
