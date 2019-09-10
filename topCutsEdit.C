void topCutsEdit(){
    char fnamed[100], fnames[100];
    sprintf(fnamed, "/u/jkubos/analiza/gitdir/Sigma_JK/outputs/anaSigmaOut_exp_all0719.root");
    sprintf(fnames, "/u/jkubos/analiza/gitdir/Sigma_JK/outputs/anaSigmaOut_sim_all0820.root");
    TFile *find = TFile::Open(fnamed, "READ");
    TFile *fins = TFile::Open(fnames, "READ");

    TH1F *hMtdd = (TH1F*)find -> Get("hMinTrackDist");
    TH1F *hMtds = (TH1F*)fins -> Get("hMinTrackDist");

    TH1F *hdVertd = (TH1F*)find -> Get("hdVert");
    TH1F *hdVerts = (TH1F*)fins -> Get("hdVert");

    hMtds -> Scale(15);
    hdVerts -> Scale(15);

    hMtdd -> SetLineColor(1);
    hMtdd -> GetXaxis() -> SetTitleSize(.05);
    hMtdd -> GetXaxis() -> SetLabelSize(.05);
    hMtdd -> GetYaxis() -> SetTitleSize(.05);
    hMtdd -> GetYaxis() -> SetLabelSize(.05);
    hMtdd -> GetXaxis() -> SetRangeUser(-5,35);
    hMtdd -> GetXaxis() -> SetTitle("MTD_L [mm]");
    hMtdd -> GetYaxis() -> SetTitle("# [a.u.]");

    hMtds -> SetLineColor(2);
    hMtds -> SetLineStyle(9);
    hMtds -> SetLineWidth(2);
    hMtds -> GetXaxis() -> SetTitleSize(.05);
    hMtds -> GetXaxis() -> SetLabelSize(.05);
    hMtds -> GetYaxis() -> SetTitleSize(.05);
    hMtds -> GetYaxis() -> SetLabelSize(.05);
    hMtds -> GetXaxis() -> SetRangeUser(-5,35);
    hMtds -> GetXaxis() -> SetTitle("MTD_L [mm]");
    hMtds -> GetYaxis() -> SetTitle("# [a.u.]");

    hdVertd -> SetLineColor(1);
    hdVertd -> GetXaxis() -> SetTitleSize(.05);
    hdVertd -> GetXaxis() -> SetLabelSize(.05);
    hdVertd -> GetYaxis() -> SetTitleSize(.05);
    hdVertd -> GetYaxis() -> SetLabelSize(.05);
    hdVertd -> GetXaxis() -> SetRangeUser(-5,75);
    hdVertd -> GetXaxis() -> SetTitle("dVert [mm]");
    hdVertd -> GetYaxis() -> SetTitle("# [a.u.]");

    hdVerts -> SetLineColor(2);
    hdVerts -> SetLineStyle(9);
    hdVerts -> SetLineWidth(2);
    hdVerts -> GetXaxis() -> SetTitleSize(.05);
    hdVerts -> GetXaxis() -> SetLabelSize(.05);
    hdVerts -> GetYaxis() -> SetTitleSize(.05);
    hdVerts -> GetYaxis() -> SetLabelSize(.05);
    hdVerts -> GetXaxis() -> SetRangeUser(-5,75);
    hdVerts -> GetXaxis() -> SetTitle("dVert [mm]");
    hdVerts -> GetYaxis() -> SetTitle("# [a.u.]");

    TGaxis::SetMaxDigits(3);

    TLine *lineMtd = new TLine(14,0,14,3700000);
    lineMtd->SetLineWidth(2);
    lineMtd->SetLineColor(1);
    TArrow *arrowMtd = new TArrow(14, 3300000, 10, 3300000, 0.03, "<--------");
    arrowMtd->SetLineWidth(2);
    arrowMtd->SetLineColor(1);
    arrowMtd->SetAngle(40);
    
    TLine *lineVert = new TLine(30,0,30,6500000);
    lineVert->SetLineWidth(2);
    lineVert->SetLineColor(1);
    TArrow *arrowVert = new TArrow(30, 5500000, 35, 5500000, 0.03, "-------->");
    arrowVert->SetLineWidth(2);
    arrowVert->SetLineColor(1);
    arrowVert->SetAngle(40);

    TPaveText *ptMtd = new TPaveText(0.55,0.7,0.6,0.8,"brNDC");
    ptMtd->SetBorderSize(0);
    ptMtd->SetFillColor(0);
    ptMtd->SetTextSize(0.04);
    ptMtd->SetTextColor(1);
    ptMtd->AddText("14 mm");
    
    TPaveText *ptVert = new TPaveText(0.55,0.7,0.6,0.8,"brNDC");
    ptVert->SetBorderSize(0);
    ptVert->SetFillColor(0);
    ptVert->SetTextSize(0.04);
    ptVert->SetTextColor(1);
    ptVert->AddText("30 mm");

    TPaveText *pt1 = new TPaveText(0.55,0.55,0.85,0.65,"brNDC");
    pt1->SetBorderSize(0);
    pt1->SetFillColor(0);
    pt1->SetTextSize(0.035);
    pt1->SetTextColor(1);
    pt1 -> AddText("data");
    pt1 -> AddText("p(3.5GeV)+p#rightarrow#Sigma(1385)^{+}+X");

    TPaveText *pt2 = new TPaveText(0.55,0.4,0.85,0.55,"brNDC");
    pt2->SetBorderSize(0);
    pt2->SetFillColor(0);
    pt2->SetTextSize(0.035);
    pt2->SetTextColor(2);
    pt2 -> AddText("simulations");
    pt2 -> AddText("p(3.5GeV)+p #rightarrow #Sigma(1385)^{+} + p + K^{0}");
    pt2 -> AddText("p(3.5GeV)+p #rightarrow #Sigma(1385)^{+} + n + K^{+}");
    


    TCanvas *cMtd = new TCanvas("cMtd", "Minimal Tracks Distance p-#pi^{-}");
    cMtd -> cd();
    hMtds -> Draw();    
    hMtdd -> Draw("same");
    lineMtd -> Draw("same");
    arrowMtd -> Draw("same");
    ptMtd -> Draw("same");
    pt1 -> Draw("same");
    pt2 -> Draw("same");
    
    TCanvas *cVert = new TCanvas("cVert", "Distance between #Sigma and #Lambda decay verticies");
    cVert -> cd();
    hdVertd -> Draw();
    hdVerts -> Draw("same");    
    lineVert -> Draw("same");
    arrowVert -> Draw("same");
    ptVert -> Draw("same");
    pt1 -> Draw("same");
    pt2 -> Draw("same");

}
