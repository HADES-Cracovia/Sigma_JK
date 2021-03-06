#include <TTree.h>


void sideBand(){
    double lm1 = 1060;
    double lm2 = 1300;
    double sm1 = 1200;
    double sm2 = 2000;
    double nbinsm = 80;
    
    TH1F *hLm_mtdL_dvert = new TH1F("hLm_mtdL_dvert", "Lambda mass after MTD p-#pi^{-} and dVert", nbinsm, lm1, lm2);
    TH1F *hSm_mtdL_dvert = new TH1F("hSm_mtdL_dvert", "Sigma mass after MTD #Lambda-#pi^{+}, dVert", nbinsm, sm1, sm2);
    TH1F *hSm_mtdL_dvert_Lm = new TH1F("hSm_mtdL_dvert_Lm", "Sigma mass after MTD #Lambda-#pi^{+}, dVert and Lmass", nbinsm, sm1, sm2);

//    hLm_mtdL_dvert -> Sumw2();
    hSm_mtdL_dvert -> Sumw2();
    
    //side band
    TH1F *hLm_sb_bg = new TH1F("hLm_sb_bg", "hLm_sb_bg", nbinsm, lm1, lm2);
    TH1F *hLm_sb_sig = new TH1F("hLm_sb_sig", "hLm_sb_sig", nbinsm, lm1, lm2);
    TH1F *hSm_sb_bg = new TH1F("hSm_sb_bg", "hSm_sb_bg", nbinsm, sm1, sm2);
    TH1F *hSm_sb_sig = new TH1F("hSm_sb_sig", "hSm_sb_sig", nbinsm, sm1, sm2);

//    hLm_sb_bg -> Sumw2();
//    hLm_sb_sig -> Sumw2();
    hSm_sb_bg -> Sumw2();
    hSm_sb_sig -> Sumw2();
    
//    double sbl1 = 1085;
    double sbl1 = 1095;
    //  double sbl2 = 1100;
    double sbl2 = 1105;
//    double sbp1 = 1130;
    double sbp1 = 1125;
//    double sbp2 = 1150;
    double sbp2 = 1140;
    //def top cuts
    float cmtd = 14; //10; //mm
    float cdvert = 30; //mm distance between verticies
    float cmm1 = 1350; //miss mass ABC, MeV
    float cmm2 = 1250; //inv mass AC, MeV                 
    float cLm1 = 1105; //MeV
    float cLm2 = 1125; //MeV
    
    TChain tree("T");
    for(int i = 106; i < 128; i++){
	char inName[100];
	printf("apr07_day_%d.root\n",i);
	sprintf(inName, "/u/jkubos/analiza/gitdir/Sigma_JK/out_packed_3/apr07_day_%d.root", i);
//	sprintf(inName, "/u/jkubos/analiza/gitdir/Sigma_JK/out_packed_3/apr07_day_%d.root", 112);
        tree.Add(inName);
    }
    
    Long_t event;
    Long_t nentries = tree.GetEntries();
    Float_t Lm, Sm, primVertX, LDecVertX, primVertY, LDecVertY, primVertZ, LDecVertZ,
	mtd, mmabc;
    tree.SetBranchAddress("fLambda_M", &Lm);
    tree.SetBranchAddress("fSigma_M", &Sm);
    tree.SetBranchAddress("fPrimaryVertexX", &primVertX);
    tree.SetBranchAddress("fPrimaryVertexY", &primVertY);
    tree.SetBranchAddress("fPrimaryVertexZ", &primVertZ);
    tree.SetBranchAddress("fLambdaDecayX", &LDecVertX);
    tree.SetBranchAddress("fLambdaDecayY", &LDecVertY);
    tree.SetBranchAddress("fLambdaDecayZ", &LDecVertZ);
    tree.SetBranchAddress("fMinTrackDist", &mtd);
    tree.SetBranchAddress("fMMABC", &mmabc);

    int eventNo = -1;
    TH1I *heventNo = new TH1I();

    heventNo -> SetName("heventNo");
    for (event = 0; event < nentries; event++){
	if(!(event%10000)){
	    printf("Event: %ld\r",event);
	}
	tree.GetEntry(event);

	float dVert = TMath::Sqrt( (LDecVertX - primVertX)*(LDecVertX - primVertX) + (LDecVertY - primVertY)*(LDecVertY - primVertY) + (LDecVertZ - primVertZ)*(LDecVertZ - primVertZ) );
	
	if(mmabc < cmm1)
	    continue;
	

	if((mtd < cmtd) && (dVert > cdvert)){
		hLm_mtdL_dvert -> Fill(Lm);
		hSm_mtdL_dvert -> Fill(Sm);
		if((Lm > sbl1 && Lm < sbl2) || (Lm > sbp1 && Lm < sbp2)){
		    hLm_sb_bg -> Fill(Lm);
		    hSm_sb_bg -> Fill(Sm);
		}else if(Lm > sbl2 && Lm < sbp1){
		    hLm_sb_sig -> Fill(Lm);
		    hSm_sb_sig -> Fill(Sm);
		}
		if(Lm > cLm1 && Lm < cLm2)
		    hSm_mtdL_dvert_Lm -> Fill(Sm); //should be the same as hSm_sb_sig
	}
	eventNo = event;
	heventNo -> Fill(eventNo);
    }//end event loop
    printf("\n");

    ///
    int a = 1090;
    int b = 1160;
    TF1 *ftmp, *fsigtmp, *fbgtmp;
    TH1F *htmp;
    double partmp[15];

    double par0 = hLm_mtdL_dvert -> GetMaximum();
    
    ftmp = new TF1("ftmp", "[0]*TMath::Voigt(x+[1],[2],[3]) + pol6(4)", a, b);
    fsigtmp = new TF1("fsigtmp", "[0]*TMath::Voigt(x+[1],[2],[3])", a, b);
    fbgtmp = new TF1("fbgtmp", "pol6(0)", a, b);
    ftmp -> SetParameters(
	//1.46482e+04,
	par0, -1.11489e+03,1.35622e+00,1.61722e+00,
	-4.21492e+06,3.65727e+03,3.32732e+00,6.71156e-05,-2.55953e-06,-2.31855e-09,2.01121e-12
	);
    htmp = (TH1F*)hLm_mtdL_dvert -> Clone("htmp");
    htmp -> Fit("ftmp", "0", "", a, b);
    htmp -> Fit("ftmp", "0", "", a, b);
    ftmp -> GetParameters(partmp);
    fsigtmp -> SetParameters(partmp);
    fbgtmp -> SetParameters(&partmp[4]);

    ftmp -> SetLineColor(kBlack);
    fsigtmp -> SetLineColor(kGreen);
    fbgtmp -> SetLineColor(kRed);

//    TH1F *hLsig = (TH1F*)hLm_mtdL_dvert -> Clone("hLsig");
    TH1F *hLsig = (TH1F*)hLm_sb_sig -> Clone("hLsig");
    hLsig -> Sumw2(); 
    hLsig -> Add(fbgtmp, -1);
    hLsig -> SetLineColor(kGreen);
    hLsig -> SetMarkerColor(kGreen);
    hLsig -> SetMarkerStyle(20);

//    TH1F *hLbg = (TH1F*)hLm_mtdL_dvert -> Clone("hLbg");
    TH1F *hLbg = (TH1F*)hLm_sb_sig -> Clone("hLbg");
    hLbg -> Sumw2();
    hLbg -> Add(hLsig, -1);
    hLbg -> SetLineColor(kRed);
    hLbg -> SetMarkerColor(kRed);
    hLbg -> SetMarkerStyle(20);

//
    double sbl = fbgtmp -> Integral(sbl1, sbl2); //A
    double sbp = fbgtmp -> Integral(sbp1, sbp2); //B
    double sb = sbl; //A+B
    double sg = fsigtmp -> Integral(sbl2, sbp1); //D
    double bg = fbgtmp -> Integral(sbl2, sbp1); //C

    double x = bg/sb;
//
    int sbl1_nb = hLm_sb_bg -> GetXaxis() -> FindBin(sbl1);
    int sbl2_nb = hLm_sb_bg -> GetXaxis() -> FindBin(sbl2);
    int sbp1_nb = hLm_sb_bg -> GetXaxis() -> FindBin(sbp1);
    int sbp2_nb = hLm_sb_bg -> GetXaxis() -> FindBin(sbp2);
    int sbg1_nb = hLm_sb_sig -> GetXaxis() -> FindBin(sbl2);
    int sbg2_nb = hLm_sb_sig -> GetXaxis() -> FindBin(sbp1);
    int sig1_nb = hLsig -> GetXaxis() -> FindBin(sbl2);
    int sig2_nb = hLsig -> GetXaxis() -> FindBin(sbp1);
    int bg1_nb = hLbg -> GetXaxis() -> FindBin(sbl2);
    int bg2_nb = hLbg -> GetXaxis() -> FindBin(sbp1);

    cout << "sbl1_nb, sbl2_nb: " << sbl1_nb << " " << sbl2_nb << " sbp1_nb, sbp2_nb: " << sbp1_nb << " " << sbp2_nb<< " sbg1_nb, sbg2_nb: " << sbg1_nb << " " << sbg2_nb << " sig1_nb, sig2_nb: " << sig1_nb << " " << sig2_nb  << " bg1_nb, bg2_nb: " << bg1_nb << " " << bg2_nb << endl;
	
    double sblh = hLm_sb_bg -> Integral(); //A
//    double sbph = hLm_sb_bg -> Integral(); //B
    double sbh = sblh; //A+B
    double sgh = hLsig -> Integral(sig1_nb, sig2_nb); //D
    double bgh = hLbg -> Integral(bg1_nb, bg2_nb); //C
    double sbgh = hLm_sb_sig -> Integral(sbg1_nb, sbg2_nb); //C+D
    
    double xh = bgh/sbh;

    hLsig -> GetXaxis() -> SetRangeUser(sbl2, sbp1);
    hLbg -> GetXaxis() -> SetRangeUser(sbl2, sbp1);
    
    cout << "integral of fitted functions:" << endl <<  " A+B=" << sb << " C=" << bg << " D=" << sg << " x=bg/sb=" << x << endl;
    cout << "integral of histograms:" << endl <<  " A+B=" << sbh << " C=" << bgh << " D=" << sgh << " C+D=" << sbgh << " xh=bgh/sbh=" << xh << endl;  

    ofstream write("countsSigma.txt", ios_base::app);
    write <<  "\n integral of fitted functions:" << endl <<  " A+B=" << sb << " C=" << bg << " D=" << sg <<  " x=bg/sb=" << x << endl;
    write << "integral of histograms:" << endl <<  " A+B=" << sbh << " C=" << bgh << " D=" << sgh << " C+D=" << sbgh << " xh=bgh/sbh=" << xh << endl;  
    
    TLegend *l1 = new TLegend(.6,.4,.8,.65);
    l1 -> SetTextSize(.035);
    l1 -> SetBorderSize(0);
    l1 -> SetFillColor(0);
    l1 -> AddEntry(hLm_mtdL_dvert, " p#pi^{-}, Cuts: MTD_L,dVert");
    l1 -> AddEntry(hLm_sb_sig, "p#pi^{-} #Lambda peak region");
    l1 -> AddEntry(hLm_sb_bg, "p#pi^{-} Sideband region");
    l1 -> AddEntry(hLsig, "#Lambda signal");
    l1 -> AddEntry(hLbg, "BG");

    hLm_sb_sig -> SetLineColor(8);
    hLm_sb_sig -> SetFillStyle(3004);
    hLm_sb_sig -> SetFillColor(kGreen);
    hLm_sb_bg -> SetLineColor(kRed);
    hLm_sb_bg -> SetFillStyle(3004);
    hLm_sb_bg -> SetFillColor(kRed);
    
    TCanvas *c1 = new TCanvas("c1","c1", 1200, 800);
    c1 -> cd();
    hLm_mtdL_dvert -> GetXaxis() -> SetTitle("M_{p#pi^{-}}");
    hLm_mtdL_dvert -> GetYaxis() -> SetTitle("a.u. [#]");
    hLm_mtdL_dvert -> GetXaxis() -> SetTitleSize(.05);
    hLm_mtdL_dvert -> GetXaxis() -> SetLabelSize(.05);

    hLm_mtdL_dvert -> Draw();
    hLm_sb_sig -> Draw("same");
    hLm_sb_bg -> Draw("same");
    hLsig -> Draw("same p");
    hLbg -> Draw("same p");
    l1 -> Draw("same");
    
//    htmp -> Draw("same");
    ftmp -> Draw("same");
    fsigtmp -> Draw("same");
    fbgtmp -> Draw("same");

    ///
    TH1F *hsbgnorm = hSm_sb_bg -> Clone("hsbgnorm");
    hsbgnorm -> Sumw2();
    hsbgnorm -> Scale(x);
    hsbgnorm -> SetLineColor(kRed);

    TH1F *hsbgnormh = hSm_sb_bg -> Clone("hsbgnormh");
    hsbgnormh -> Sumw2();
    hsbgnormh -> Scale(xh);
    hsbgnormh -> SetLineColor(4);
    hsbgnormh -> SetMarkerColor(4);
    hsbgnormh -> SetMarkerStyle(20);

    ///
    hSm_sb_sig -> SetLineColor(8);

    TH1F *hsignal1 = hSm_sb_sig -> Clone("hsignal");
    hsignal1 -> Sumw2();
    hsignal1 -> Add(hsbgnorm, -1);
    hsignal1 -> SetLineColor(kMagenta);

    TH1F *hsignal1h = hSm_sb_sig -> Clone("hsignalh");
    hsignal1h -> Sumw2();
    hsignal1h -> Add(hsbgnormh, -1);
    hsignal1h -> SetLineColor(kMagenta);
    hsignal1h -> SetMarkerColor(kMagenta);
    hsignal1h -> SetMarkerStyle(20);

    //fit
    double par1 = hsignal1h -> GetMaximum();
    
    int fitS1 = 1255;
    int fitS2 = 1560;
    
    TF1 * fhistSmtdLDvertLm = new TF1("fhistSmtdLDvertLm", "[0]*[1]*[1]/( (x*x-[2]*[2])*(x*x-[2]*[2]) + (x*x*x*x*[1]*[1]/([2]*[2])) ) + pol6(3)", fitS1, fitS2);
    fhistSmtdLDvertLm -> SetParameters(
	par1,40,1376.04,-747286,-3544.78,9.34086,-0.00679285,1.57678e-06,6.18276e-16
	);
    fhistSmtdLDvertLm -> SetParLimits(0,0,1000000);
    fhistSmtdLDvertLm -> SetParLimits(1,30,45);
    fhistSmtdLDvertLm -> SetParLimits(2,1360,1400);

    hsignal1h -> Fit("fhistSmtdLDvertLm", "0", "", fitS1, fitS2);
    hsignal1h -> Fit("fhistSmtdLDvertLm", "0", "", fitS1, fitS2);

    TF1 * fsigSmtdLDvertLm = new TF1("fsigSmtdLDvertLm", "[0]*[1]*[1]/( (x*x-[2]*[2])*(x*x-[2]*[2]) + (x*x*x*x*[1]*[1]/([2]*[2])) )", fitS1, fitS2);
    TF1 * fbgSmtdLDvertLm = new TF1("fbgSmtdLDvertLm", "pol6(0)", fitS1, fitS2);
    double parSmtdLDvertLm[12];
    fhistSmtdLDvertLm -> GetParameters(parSmtdLDvertLm);
    fsigSmtdLDvertLm -> SetParameters(parSmtdLDvertLm);
    fbgSmtdLDvertLm -> SetParameters(&parSmtdLDvertLm[3]);

    fsigSmtdLDvertLm -> SetLineColor(kBlack);
    fsigSmtdLDvertLm -> SetMarkerColor(kBlack);
    fsigSmtdLDvertLm -> SetMarkerStyle(20);
    fsigSmtdLDvertLm -> SetMarkerSize(.5);
    fbgSmtdLDvertLm -> SetLineColor(kBlue);
    fbgSmtdLDvertLm -> SetMarkerColor(kBlue);
    fbgSmtdLDvertLm -> SetMarkerStyle(20);
    fbgSmtdLDvertLm -> SetMarkerSize(.5);
    
    int sigmaAll1_nb = hSm_sb_sig -> GetMinimumBin();
    int sigmaAll2_nb = hSm_sb_sig -> GetMaximumBin();
    int sigmaSB1_nb = hSm_sb_bg -> GetMinimumBin();
    int sigmaSB2_nb = hSm_sb_bg -> GetMaximumBin();
    int sigmaSBnorm1_nb = hsbgnormh -> GetMinimumBin();
    int sigmaSBnorm2_nb = hsbgnormh -> GetMaximumBin();
    int sigma1_nb = hsignal1h -> GetMinimumBin();
    int sigma2_nb = hsignal1h -> GetMaximumBin();
    
    double intSigma = (hsignal1->GetBinWidth(0))*(hsignal1 -> Integral());
    cout << "Sigma=" << intSigma << endl;
    write << "Sigma=" << intSigma << endl;

    
    double intSigmahAll = hSm_sb_sig -> Integral();
    double intSigmahSB = hSm_sb_bg -> Integral();
    double intSigmahSBnorm = hsbgnormh -> Integral();
    double intSigmah = hsignal1h -> Integral();

    cout << "SigmaAll=" << intSigmahAll << " SigmaSB=" << intSigmahSB << " SigmaSBnorm=" << intSigmahSBnorm << " Sigmah=" << intSigmah << endl;
    write << "SigmaAll=" << intSigmahAll << " SigmaSB=" << intSigmahSB << " SigmaSBnorm=" << intSigmahSBnorm << " Sigmah=" << intSigmah << endl;
    write.close();

    hSm_sb_bg -> SetLineColor(kRed);
//    hSm_sb_bg -> SetLineStyle(9);

    TLegend *l2 = new TLegend(.55,.4,.7,.65);
    l2 -> SetTextSize(.035);
    l2 -> SetBorderSize(0);
    l2 -> SetFillColor(0);
    l2 -> AddEntry(hSm_sb_sig, "1. p#pi^{-} from #Lambda peak region + #pi^{+}");
    l2 -> AddEntry(hSm_sb_bg, "2. p#pi^{-} from Sideband region + #pi^{+}");
//    l2 -> AddEntry(hsbgnorm, "3a. p#pi^{-} from Sideband region + #pi^{+} (norm, fit)");
//    l2 -> AddEntry(hsignal1, "1 - 3a");
    l2 -> AddEntry(hsbgnormh, "3. p#pi^{-} from Sideband region + #pi^{+} (norm)");
    l2 -> AddEntry(hsignal1h, "1 - 3");

    TPaveText *pt = new TPaveText(.55,.7,.6,.85, "NDC");
    pt -> SetTextSize(.033);
    pt -> SetFillColor(0);
    pt -> SetBorderSize(0);
    char nCD[32], nAB[64], nC[32], nD[32],  nstsb[32];
    sprintf(nCD, "N Signal region = %.0f", sbgh);
    sprintf(nC, "N BG in Signal region = %.0f", bgh);
    sprintf(nD, "N S in Signal region = %.0f", sgh);
    sprintf(nAB, "N Sideband = %.0f", sblh);
//    sprintf(stsb, "StSB = %.3f", StSB);
    pt -> AddText(nCD);
    pt -> AddText(nAB);
    pt -> AddText(nC);
    pt -> AddText(nD);
//    pt -> AddText(stsb);

    TPaveText *pt2 = new TPaveText(.55,.7,.6,.85, "NDC");
    pt2 -> SetTextSize(.033);
    pt2 -> SetFillColor(0);
    pt2 -> SetBorderSize(0);
    char nCD2[64], nAB2[64], nC2[64], nD2[64];
    sprintf(nCD2, "#Sigma from #Lambda Signal region = %d", intSigmahAll);
    sprintf(nAB2, "BG from Sideband = %d", intSigmahSB);
    sprintf(nC2, "BG from Signal region = %d", intSigmahSBnorm);
    sprintf(nD2, "#Sigma from #Lambda peak = %d", intSigmah);
    pt2 -> AddText(nCD2);
    pt2 -> AddText(nAB2);
    pt2 -> AddText(nC2);
    pt2 -> AddText(nD2);
    

    TCanvas *c2 = new TCanvas("c2","c2",1200,800);
    c2 -> cd();
    hSm_sb_sig -> GetXaxis() -> SetTitle("M_{p#pi^{-}#pi^{+}}");
    hSm_sb_sig -> GetYaxis() -> SetTitle("a.u. [#]");
    hSm_sb_sig -> GetXaxis() -> SetTitleSize(.05);
    hSm_sb_sig -> GetXaxis() -> SetLabelSize(.05);
    
    hSm_sb_sig -> Draw("same");
    hSm_sb_bg -> Draw("same");
//    hsbgnorm -> Draw("same");
//    hsignal1 -> Draw("same");
    hsbgnormh -> Draw("same p");
    hsignal1h -> Draw("same p");
    fsigSmtdLDvertLm ->  Draw("same");
    fbgSmtdLDvertLm -> Draw("same");
    l2 -> Draw("same");
    pt2 -> Draw("same");
    
    c1 -> cd();
    pt -> Draw("same");
    
    TFile *fout = TFile::Open("./sb_out_3_test.root", "RECREATE");
    hLm_mtdL_dvert -> Write();
    hLm_sb_bg -> Write();
    hLm_sb_sig -> Write();
    hSm_sb_bg -> Write();
    hSm_sb_sig -> Write();
    hSm_mtdL_dvert -> Write();
    hSm_mtdL_dvert_Lm -> Write();
    hLsig -> Write();
    hLbg -> Write();
    c1 -> Write();
    hsbgnorm -> Write();
    hsignal1 -> Write();
    hsbgnormh -> Write();
    hsignal1h -> Write();
    c2 -> Write();
    
    heventNo -> Write();

    fout -> Close();
        


}
