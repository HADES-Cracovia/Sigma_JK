TH1F* bgSubtr(TH1F* h1){

  TH1F *hist = (TH1F*)h1 -> Clone();
  const char *hname = hist -> GetName();
  char namePeak[50];
  sprintf(namePeak, "%s_peak", hname);

  cout << ">>>>>>>>fit BG<<<<<<<<<" << endl;
  TF1 * fhist = new TF1("fhist", "gaus(0)+gaus(3)+pol1(6)", 1080,1140);
  fhist -> SetParameters(
	10000,1.11421e+03,5,
	10000,1.11450e+03,5,
	2.5);
	//2984560,1114.29,4.12196,
	//8.78103e+06,1.11469e+03,1.95058,
	//-4.32460e+09,7.51569e+06,-3.25342e+03,1,1);
  fhist -> SetParLimits(0, 0, 100000000);
  fhist -> SetParLimits(1, 1110, 1120);
  fhist -> SetParLimits(2, 0, 10);
  fhist -> SetParLimits(3, 0, 100000000);
  fhist -> SetParLimits(4, 1110, 1120);
  fhist -> SetParLimits(5, 0, 10);
  hist -> Fit("fhist", "", "", 1080, 1180);
  hist -> Fit("fhist", "", "", 1080, 1180);
  hist -> SetName(namePeak);
  TF1 * fsig = new TF1("fsig", "gaus(0)+gaus(3)", 1080,1140);
  TF1 * fbg = new TF1("fbg", "pol1(0)", 1080,1140);
  double par[12];
  fhist -> GetParameters(par);
  fsig -> SetParameters(par);
  fbg -> SetParameters(&par[6]);
  
  TH1F *hout = (TH1F*)hist -> Clone("histBG");
  hout -> Add(fsig, -1);

  hout -> SetLineColor(kRed);
  
  return hout;

}
