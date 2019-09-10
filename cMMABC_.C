void cMMABC_(){
    char fnamed[100], fnames[100];
    sprintf(fnamed, "/u/jkubos/analiza/gitdir/Sigma_JK/outputs/anaSigmaOut_exp_all0625.root");
    sprintf(fnames, "/u/jkubos/analiza/gitdir/Sigma_JK/outputs/anaSigmaOut_sim_all0625.root");
    TFile *find = TFile::Open(fnamed, "READ");
    TFile *fins = TFile::Open(fnames, "READ");
    
    TH1F *hMMd = (TH1F*)find -> Get("hMMABC0");
    TH1F *hMMs = (TH1F*)fins -> Get("hMMABC0");
    hMMs -> Scale(20);

    hMMd -> SetLineColor(1);
    hMMd -> GetXaxis() -> SetTitleSize(.05);
    hMMd -> GetXaxis() -> SetLabelSize(.05);
    hMMd -> GetYaxis() -> SetTitleSize(.05);
    hMMd -> GetYaxis() -> SetLabelSize(.05);
    hMMd -> GetXaxis() -> SetRangeUser(0,2500);
    hMMd -> GetXaxis() -> SetTitle("MM_{p#pi^{-}#pi^{+}} [MeV]");
    hMMd -> GetYaxis() -> SetTitle("# [a.u.]");
    TGaxis::SetMaxDigits(3);
    hMMd -> SetTitle("p#pi^{-}#pi^{+} missing mass");
    
    hMMs -> SetLineColor(2);
    hMMs -> SetLineStyle(9);
    hMMs -> SetLineWidth(2);
    hMMs -> GetXaxis() -> SetTitleSize(.05);
    hMMs -> GetXaxis() -> SetLabelSize(.05);
    hMMs -> GetYaxis() -> SetTitleSize(.05);
    hMMs -> GetYaxis() -> SetLabelSize(.05);

    TLine *line = new TLine(1350,0,1350,10000000);
    line->SetLineWidth(2);
    line->SetLineColor(1);
    TArrow *arrow = new TArrow(1350, 9000000, 1750, 9000000, 0.03, "-------->");
    arrow->SetLineWidth(2);
    arrow->SetLineColor(1);
    arrow->SetAngle(40);

    TPaveText *pt1 = new TPaveText(0.24,0.65,0.27,0.8,"brNDC");
    pt1->SetBorderSize(0);
    pt1->SetFillColor(0);
    pt1->SetTextSize(0.033);
    pt1->SetTextColor(1);
    pt1 -> AddText("data");
    pt1 -> AddText("p(3.5GeV)+p#rightarrow#Sigma(1385)^{+}+X");

    TPaveText *pt2 = new TPaveText(0.72,0.5,0.74,0.7,"brNDC");
    pt2->SetBorderSize(0);
    pt2->SetFillColor(0);
    pt2->SetTextSize(0.033);
    pt2->SetTextColor(2);
    pt2 -> AddText("simulations");
    pt2 -> AddText("p(3.5GeV)+p #rightarrow #Sigma(1385)^{+} + p + K^{0}");
    pt2 -> AddText("p(3.5GeV)+p #rightarrow #Sigma(1385)^{+} + n + K^{+}");

    TPaveText *pt3 = new TPaveText(0.55,0.83,0.6,0.83,"brNDC");
    pt3->SetBorderSize(0);
    pt3->SetFillColor(0);
    pt3->SetTextSize(0.033);
    pt3->SetTextColor(1);
    pt3 -> AddText("1350 MeV");

    TCanvas *c1 = new TCanvas("cMMABC_", "p#pi^{-}#pi^{+} missing mass");
    c1->cd();
    hMMd -> Draw();
    hMMs -> Draw("same");
    line -> Draw("same");
    arrow -> Draw("same");
    pt1 -> Draw("same");
    pt2 -> Draw("same");
    pt3 -> Draw("same");
}
