void peakBG_()
{
//=========Macro generated from canvas: Canvas_1/Canvas_1
//=========  (Wed Jul 10 11:07:50 2019) by ROOT version 6.14/06
  TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1", 1200, 800);

    char fnamed[100], fnames[100];
    sprintf(fnamed, "./outputs/anaSigmaOut_exp_all0919_3.root");
    sprintf(fnames, "./outputs/anaSigmaOut_sim_all0820.root");
    TFile *find = TFile::Open(fnamed, "READ");
    TFile *fins = TFile::Open(fnames, "READ");

    TH1F *hSm_mtdL_dvert_Lm__1 = (TH1F*)find -> Get("hSm_mtdL_dvert_Lm");
    TH1F *hSm_mtdL_dvert_Lm__2 = (TH1F*)fins -> Get("hSm_mtdL_dvert_Lm");
    hSm_mtdL_dvert_Lm__2->Scale(.08);
    hSm_mtdL_dvert_Lm__2->SetMarkerStyle(5);
  
   const int n = 105;
   double x[n], y[n];

   int abin = hSm_mtdL_dvert_Lm__1 -> GetMinimumBin();
   int bbin = hSm_mtdL_dvert_Lm__1 -> GetMaximumBin();
   double a = hSm_mtdL_dvert_Lm__1 -> GetBinCenter(abin);
   double b = hSm_mtdL_dvert_Lm__1 -> GetBinCenter(bbin);
   double a1 = 1400.;
   double b1 = 1500.;
   double min = 1250.;
   double max = 1800.;
   
   int j = abin;
   double k = 0;
   int l = 0;
   for(int i = 0; i < n; i++){
     k = hSm_mtdL_dvert_Lm__1 -> GetBinCenter(j);
     if(k < min || k > max){
       j++;
       continue;
     }
     if(k < a1 || k > b1){
       x[l] = k;
       y[l] = hSm_mtdL_dvert_Lm__1 -> GetBinContent(j);
       cout << " l x y: " << l << " " <<  x[l] << " " << y[l] << endl;
       l++;
     }
     j++;
   }

   double ymax = 1.02*(hSm_mtdL_dvert_Lm__1 -> GetBinContent(hSm_mtdL_dvert_Lm__1->GetMaximumBin()));
   TGraph *gr = new TGraph(n, x, y);
   gr -> SetTitle("inv mass #Lambda-#pi^{+}");
   gr -> GetXaxis() -> SetRangeUser(min,max);
   gr -> GetYaxis() -> SetRangeUser(0,ymax);
   gr -> GetXaxis() -> SetTitleSize(.05);
   gr -> GetXaxis() -> SetTitle("m_{p#pi^{-}#pi^{+}} [MeV]");
   gr -> GetYaxis() -> SetTitle("a.u. [#]");
   gr -> Draw("ap*");

   int fitS1 = 1255;
   int fitS2 = 1520;
   int fitSsig1 = 1325;
   int fitSsig2 = 1445;

   TF1 * fhistSmtdLDvertLm = new TF1("fhistSmtdLDvertLm", "[0]*[1]*[1]/( (x*x-[2]*[2])*(x*x-[2]*[2]) + (x*x*x*x*[1]*[1]/([2]*[2])) ) + pol6(3)", fitS1, fitS2);
   fhistSmtdLDvertLm -> SetParameters(
				      4.62946e+09,40,1376.04,-747286,-3544.78,9.34086,-0.00679285,1.57678e-06,6.18276e-16
       );
   fhistSmtdLDvertLm -> SetParLimits(0,0,10000000000);
   fhistSmtdLDvertLm -> SetParLimits(1,30,40);
   fhistSmtdLDvertLm -> SetParLimits(2,1360,1400);

   gr -> Fit("fhistSmtdLDvertLm", "0", "", fitS1, fitS2);
   gr -> Fit("fhistSmtdLDvertLm", "0", "", fitS1, fitS2);

   TF1 * fsigSmtdLDvertLm = new TF1("fsigSmtdLDvertLm", "[0]*[1]*[1]/( (x*x-[2]*[2])*(x*x-[2]*[2]) + (x*x*x*x*[1]*[1]/([2]*[2])) )", fitSsig1, fitSsig2);
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
   
   fhistSmtdLDvertLm -> Draw("same");
   fsigSmtdLDvertLm -> Draw("same");
   fbgSmtdLDvertLm -> Draw("same");
   hSm_mtdL_dvert_Lm__1->Draw("same");
   hSm_mtdL_dvert_Lm__2->Draw("same");

   TH1F *hpeakBG = (TH1F*)hSm_mtdL_dvert_Lm__1 -> Clone("hpeakBG");
   hpeakBG -> SetLineStyle(9);
   hpeakBG -> SetLineColor(8);
   hpeakBG -> SetLineWidth(2);
   hpeakBG -> Add(fhistSmtdLDvertLm, -1);
   hpeakBG -> GetXaxis() -> SetRangeUser(fitS1, fitS2);
   hpeakBG -> Draw("same");

   TLegend *leg = new TLegend(.5,.5,.85,.8);
   leg -> SetTextSize(.035);
   leg -> SetFillStyle(0);
   leg -> SetBorderSize(0);
   leg -> AddEntry(hSm_mtdL_dvert_Lm__1, "1. inv mass #Lambda-#pi^{+}");
   leg -> AddEntry(fhistSmtdLDvertLm, "2. fit BW+pol6 (w/o [1400;1500]MeV)");
   leg -> AddEntry(fsigSmtdLDvertLm, "3. fit signal (BW)");
   leg -> AddEntry(fbgSmtdLDvertLm, "4. fit BG (pol6)");
   leg -> AddEntry(hpeakBG, "1-2");
   leg -> AddEntry(hSm_mtdL_dvert_Lm__2, "simulation of 2 signal channels");
   leg -> Draw("same");


   TCanvas *Canvas_2 = new TCanvas("Canvas_2", "Canvas_2", 1200, 800);
   Canvas_2 -> cd();
   gr -> Draw("ap*");
   fhistSmtdLDvertLm -> Draw("same");
   fsigSmtdLDvertLm -> Draw("same");
   fbgSmtdLDvertLm -> Draw("same");
   hSm_mtdL_dvert_Lm__1->Draw("same");
   
   double cntSmSigSmtdLDvertLm = fsigSmtdLDvertLm -> Integral(fitSsig1,fitSsig2);
   double cntSmBGSmtdLDvertLm = fbgSmtdLDvertLm -> Integral(fitSsig1,fitSsig2);
   char textbw[64], gam[64], mSig[64], nSig[64];
   sprintf(textbw, "Breit-Wigner:\n");
   sprintf(gam, "#Gamma = %.2f MeV", parSmtdLDvertLm[1]);
   sprintf(mSig, "M = %.2f", parSmtdLDvertLm[2]);
   sprintf(nSig, "N = %.0f", cntSmSigSmtdLDvertLm);
   TPaveText *ptBW = new TPaveText(.75, .4, .8, .6, "NDC");
   ptBW -> SetFillColor(0);
   ptBW -> SetBorderSize(0);
   ptBW -> SetTextSize(0.04);
   ptBW -> AddText(textbw);
   ptBW -> AddText(gam);
   ptBW -> AddText(mSig);
   ptBW -> AddText(nSig);
   ptBW -> Draw("same");    
}
