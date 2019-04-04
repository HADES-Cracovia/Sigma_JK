#include <TChain.h>
#include "TFile.h"
#include "TH1F.h"
#include <Riostream.h>
#include <TCutG.h>
#include <math.h>  
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <fstream>
//#include "TTreeReader.h"
//#include "TTreeReaderValue.h"
#include "bgSubtr.C"

void sigmaAna(){
   //def histos
   TH1F *hLm = new TH1F("hLm", "#Lambda mass", 200, 1060, 1200);
   TH1F *hpm = new TH1F("hpm", "proton mass", 200, 935, 940);
   TH1F *hpimm = new TH1F("hpimm", "#pion^{-} mass", 200, 130, 150);
   TH1F *hpipm = new TH1F("hpipm", "#pion^{+} mass", 200, 130, 150);
   TH1F *hSm = new TH1F("hSm", "#Sigma mass", 100, 1300, 1450);

   TH1F *hLth = new TH1F("hLth", "Lambda theta", 200, 0, 4);
   TH1F *hpth = new TH1F("hpth", "proton theta", 200, 0, 4);
   TH1F *hpimth = new TH1F("hpimth", "pion- theta", 200, 0, 4);
   TH1F *hpipth = new TH1F("hpipth", "pion+ theta", 200, 0, 4);
   TH1F *hSth = new TH1F("hSth", "Sigma theta", 200, 0, 4);

   TH1F *hLp = new TH1F("hLp", "Lambda momentum", 200, 0, 4000);
   TH1F *hpp = new TH1F("hpp", "proton momentum", 200, 0, 4000);
   TH1F *hpimp = new TH1F("hpimp", "pion- momentum", 200, 0, 4000);
   TH1F *hpipp = new TH1F("hpipp", "pion+ momentum", 200, -10, 20);
   TH1F *hSp = new TH1F("hSp", "Sigma momentum", 200, 0, 4000);

   TH2I *hdEdx_Mdc = new TH2I("hdEdx_Mdc", "dEdx MDC", 500, -1500, 2500, 500, 0, 10);
   TH2I *hdEdx_Mdc_acc = new TH2I("hdEdx_Mdc_acc", "dEdx MDC acc", 500, -1500, 2500, 500, 0, 10);

   TH1F *hMinTrackDist = new TH1F("hMinTrackDist", "MinTrackDist p-#pi^{-}", 100, -15, 42);
   TH1F *hMinTrackDistLpi = new TH1F("hMinTrackDistLpi", "MinTrackDist #Lambda-#pi^{+}", 100, 10, 110);
   TH1F *hPrimVertR = new TH1F("hPrimVertR", "Primary vertex R", 100, -100, 150);
   TH1F *hLDecVertR = new TH1F("hLDecVertR", "Lambda decay vertex R", 100, -100, 150);
   TH1F *hdLVR = new TH1F("hdLVR", "Lambda decay vertex - Primary Vertex", 100, -100, 150);
   TH1F *hPrimVertZ = new TH1F("hPrimVertZ", "Primary vertex Z", 100, -100, 150);
   TH1F *hLDecVertZ = new TH1F("hLDecVertZ", "Lambda decay vertex Z", 100, -100, 150);

   //mtd L cut
   TH1F *hLm_mtdL = new TH1F("hLm_mtdL", "#Lambda mass after MTD p-#pi^{-}", 200, 1060, 1200);
   TH1F *hpipm_mtdL = new TH1F("hpipm_mtdL", "#pion^{+} mass after MTD p-#pi^{-}", 200, 130, 150);
   TH1F *hSm_mtdL = new TH1F("hSm_mtdL", "#Sigma mass after MTD p-#pi^{-}", 100, 1300, 1450);

   //
   //mtd L & VertLz cut
   TH1F *hLm_mtdL_vertLz = new TH1F("hLm_mtdL_vertLz", "Lambda mass after MTD p-#pi^{-} and LDecVertZ", 200, 1060, 1200);
   TH1F *hSm_mtdL_vertLz = new TH1F("hSm_mtdL_vertLz", "Sigma mass after MTD #Lambda-#pi^{+} and LDecVertZ", 100, 1300, 1450);

   //mtd L & VertLz cut & dLV cut
   TH1F *hLm_mtdL_vertLz_dLV = new TH1F("hLm_mtdL_dLV", "Lambda mass after MTD p-#pi^{-}, LDecVertZ and (LDecVertR-PrimVertR)", 200, 1060, 1200);
   TH1F *hpipm_mtdL_vertLz_dLV = new TH1F("hpipm_mtdL_dLV", "pion+ mass after MTD p-#pi^{-}, LDecVertZ and (LDecVertR-PrimVertR)", 200, 130, 150);
   TH1F *hSm_mtdL_vertLz_dLV = new TH1F("hSm_mtdL_dLV", "Sigma mass after MTD #Lambda-#pi^{+}, LDecVertZ and (LDecVertR-PrimVertR)", 100, 1300, 1450);

   //mtd L & VertLz cut & dLV cut & Lm
   TH1F *hLm_mtdL_vertLz_dLV_Lm = new TH1F("hLm_mtdL_vertLz_dLV_Lm", "Lambda mass after MTD p-#pi^{-}, LDecVertZ, (LDecVertR-PrimVertR) and Lmass", 200, 1060, 1200);
   TH1F *hSm_mtdL_vertLz_dLV_Lm = new TH1F("hSm_mtdL_vertLz_dLV_Lm", "Sigma mass after MTD #Lambda-#pi^{+}, LDecVertZ, (LDecVertR-PrimVertR) and Lmass", 100, 1300, 1450);
   //
   //mtd L & dLV cut
   TH1F *hLm_mtdL_dLV = new TH1F("hLm_mtdL_dLV", "Lambda mass after MTD p-#pi^{-} and (LDecVertR-PrimVertR)", 200, 1060, 1200);
   TH1F *hpipm_mtdL_dLV = new TH1F("hpipm_mtdL_dLV", "pion+ mass after MTD p-#pi^{-} and (LDecVertR-PrimVertR)", 200, 130, 150)
;   TH1F *hSm_mtdL_dLV = new TH1F("hSm_mtdL_dLV", "Sigma mass after MTD #Lambda-#pi^{+} and (LDecVertR-PrimVertR)", 100, 1300, 1450);

   //mtd L & dLV cut & Lm
   TH1F *hLm_mtdL_dLV_Lm = new TH1F("hLm_mtdL_dLV_Lm", "Lambda mass after MTD p-#pi^{-}, (LDecVertR-PrimVertR) and Lmass", 200, 1060, 1200);
   TH1F *hSm_mtdL_dLV_Lm = new TH1F("hSm_mtdL_dLV_Lm", "Sigma mass after MTD #Lambda-#pi^{+}, (LDecVertR-PrimVertR) and Lmass", 50, 1300, 1450);
   //
   //dLV cut
   TH1F *hLm_dLV = new TH1F("hLm_dLV", "Lambda mass after (LDecVertR-PrimVertR)", 200,1060, 1200);
   TH1F *hpipm_dLV = new TH1F("hpipm_dLV", "pion+ mass after (LDecVertR-PrimVertR)", 200, 130, 150);
   TH1F *hSm_dLV = new TH1F("hSm_dLV", "Sigma mass after (LDecVertR-PrimVertR)", 100, 1300, 1450);

   //dLV cut & Lm
   TH1F *hLm_dLV_Lm = new TH1F("hLm_dLV_Lm", "Lambda mass after (LDecVertR-PrimVertR) and Lmass", 200, 1060, 1200);
   TH1F *hSm_dLV_Lm = new TH1F("hSm_dLV_Lm", "Sigma mass after (LDecVertR-PrimVertR) and Lmass", 100, 1300, 1450);

   //S: all L cuts &:
   //mtdS
   TH1F *hSm_mtdS = new TH1F("hSm_mtdS", "#Sigma mass after all cuts on #Lambda and MTD #Lambda-#pi^{+}", 50, 1300, 1450);
   //mtdS & primVertZ
   TH1F *hSm_mtdS_pvz = new TH1F("hSm_mtdS_pvz", "#Sigma mass after all cuts on #Lambda, MTD #Lambda-#pi^{+} and PrimVertZ", 50, 1300, 1450);
   
   //histos BG
   TH1F *hLmRecBG = new TH1F("hLmRecBG", "hLmRecBG", 200, 1060, 1200);
   TH1F *hLmRecBGMtd = new TH1F("hLmRecBGMtd", "hLmRecBGMtd", 200, 1060, 1200);
   TH1F *hLmRecBGMtdVertLz = new TH1F("hLmRecBGMtdVertLz", "hLmRecBGMtdVertLz", 200, 1060, 1200);
   TH1F *hLmRecBGMtdDlv = new TH1F("hLmRecBGMtdDlv", "hLmRecBGMtdDlv", 200, 1060, 1200);
   TH1F *hLmRecBGDlv = new TH1F("hLmRecBGDlv", "hLmRecBGDlv", 200, 1060, 1200);
//   int cntLmRec, cntLmRecBG, cntLmRecMtd, cntLmRecMtdBG, cntLmRecMtdVertLz, cntLmRecMtdVertLzBG, cntLmRecMtdDlv, cntLmRecMtdDlvBG, cntLmRecDlv, cntLmRecDlvBG;
   //canvases with S&BG
   TCanvas *cLsbg = new TCanvas("cLsbg", "cLsbg", 1200, 800);
   TCanvas *cSsbg = new TCanvas("cSsbg", "cSsbg", 1200, 800);
   cLsbg -> Divide(2,2);
   cSsbg -> Divide(2,2);
   
   //miss mass spectra
   TH1F *hMMABC = new TH1F("hMMABC", "hMMABC", 100, -2300, 2350);
   TH1F *hMMABCpip = new TH1F("hMMABCpip", "hMMABCpip", 100, -3600, 2400);
   TH1F *hMMABCpim = new TH1F("hMMABCpim", "hMMABCpim", 100, -3600, 2400);
   TH1F *hMMAB = new TH1F("hMMAB", "hMMAB", 100, -2300, 2350);
   TH1F *hMMBC = new TH1F("hMMBC", "hMMBC", 100, -2300, 2350);
   //inv mass spectra
   TH1F *hinvMAC = new TH1F("hinvMAC", "hinvMAC", 100, 1055, 2000);
   TH1F *hinvMAB = new TH1F("hinvMAB", "hinvMAB", 100, 1066, 1212);

   TH2F *hDN = new TH2F("hDN", "hDN", 100, 1066, 1212, 100, -3600, 2400);


   //read graph cuts
   char fcutp[100], fcutpip[100], fcutpim[100];
   char cutp[50] = "Mdc_dEdx_P_cut_mod_ChiiV1"; 
   char cutpip[50] = "Mdc_dEdx_PiP_cut_PID_mod_ChiiV2";
   char cutpim[50] = "Mdc_dEdx_PiM_cut_PID_mod_ChiiV2";
   //char cutp[50] = "CUTG_14"; 
   //char cutpip[50] = "CUTG_8";
   //char cutpim[50] = "CUTG_9";
   sprintf(fcutp, "%s.root", cutp);
   sprintf(fcutpip, "%s.root", cutpip);
   sprintf(fcutpim, "%s.root", cutpim);
   cout << "cuts: " << fcutp << " " << fcutpip << " " << fcutpim << endl;
   TFile *fCutp = TFile::Open(fcutp);
   TFile *fCutpip = TFile::Open(fcutpip);
   TFile *fCutpim = TFile::Open(fcutpim);
   TCutG *gCutp, *gCutpip, *gCutpim;
   gCutp = (TCutG*)fCutp -> Get(cutp);
   gCutpip = (TCutG*)fCutpip -> Get(cutpip);
   gCutpim = (TCutG*)fCutpim -> Get(cutpim);

   fCutp -> Close();
   fCutpip -> Close();
   fCutpim -> Close();
   
   /*   //myCuts
   TFile *fCuts = TFile::Open("myCuts.root");   
   int partN = 15;
   TCutG *gCut[partN];
   int k = 0;
   int isCut[15];
   for(int pidNo = 0; pidNo < partN; pidNo++){
     char cutName[50];
     sprintf(cutName, "CUTG_%d", pidNo);
     if(!(TCutG*)fCuts -> Get(cutName)){
       sprintf(cutName, "CUTG_%d", 0);
       continue;
     }
     gCut[pidNo] = (TCutG*)fCuts -> Get(cutName);
     isCut[pidNo] = 1;
     cout << "cutName: " << cutName << endl;
     k++;
   }
   fCuts -> Close();
   *///end myCuts

   //def top cuts
   float cmtd = 10; //mm
   float cvertLz1 = -20; //mm
   float cvertLz2 = 300; //mm
   float csVL1 = 20; //mm
   float csVL2 = 80; //mm
   float cLm1 = 1105; //MeV
   float cLm2 = 1125; //MeV
   float cmtdLpi = 40; //mm
   float cpvz1 = -55; //mm
   float cpvz2 = -30; //mm
   float cmm1 = 1050; //miss mass ABC, MeV
   float cmm2 = 1150; //miss mass ABC+pip, MeV
   
   /* //variables for the MTD L scan
   int mtdLi = 0;
   TH1F *hmtdLscan[15];
   TCanvas *cMtdScan = new TCanvas("cMtdScan", "cMtdScan", 1200, 800);
   cMtdScan -> Divide(5,3);
   char nameMtdLi[30];
   double mtdVal[15], sgnVal[15];
   double sg, bg;
   for(int i = 0; i < 15; i++){
       mtdLi = (i+1)*2;
       sprintf(nameMtdLi, "hmtdLscan_%d", mtdLi);
       hmtdLscan[i] = new TH1F("hmtdLscan", nameMtdLi, 200, 1060, 1200);
   }
   //variables for the dLV scan
   int dLVi = 0;
   TH1F *hdLVscan[15];
   TCanvas *cdLVScan = new TCanvas("cdLVScan", "cdLVScan", 1200, 800);
   cdLVScan -> Divide(5,3);
   char namedLVi[30];
   double dLVVal[15], sgnVal2[15];
   double sg2, bg2;
   for(int i = 0; i < 15; i++){
       dLVi = (i+1)*5;
       sprintf(namedLVi, "hdLVscan_%d", dLVi);
       hdLVscan[i] = new TH1F("hdLVscan", namedLVi, 200, 1060, 1200);
   }
   //variables for the MtdL(dLV) scan
   int dLVi_mtdL = 0;
   int mtdLi_dLV = 0;
   TCanvas *cmtdLdLV[15];
   TH1F *hdLVmtdLscan[15][15];
   TGraph2D *gdLVmtdL = new TGraph2D();
   char namedLVi_mtdL[30];
   char nameCanv[30];
   double sgnVal3[15][15];
   double sg3, bg3;
   for(int i = 0; i < 15; i++){
       sprintf(nameCanv, "cmtdLdLV_%d", i);
       cmtdLdLV[i] = new TCanvas(nameCanv, nameCanv, 1200, 800);
       cmtdLdLV[i] -> Divide(5,3);
       mtdLi_dLV = (i+1)*2;
       for(int j = 0; j < 15; j++){
	   dLVi_mtdL = (j+1)*5;
	   sprintf(namedLVi_mtdL, "hdLV_%d_mtdL_%d_scan", dLVi_mtdL, mtdLi_dLV);
	   hdLVmtdLscan[i][j] = new TH1F("hdLVmtdLscan", namedLVi_mtdL, 200, 1060, 1200);
       }
   }
   */
   //variables for the MTD S scan                                                                                 
   int mtdLpii = 0;
   TH1F *hmtdLpiscan[12];
   TCanvas *cMtdLpiScan = new TCanvas("cMtdLpiScan", "cMtdLpiScan", 1200, 800);
   cMtdLpiScan -> Divide(4,3);
   char nameMtdLpii[30];
   double mtdLpiVal[12], sgnVal4[12];
   double sg4, bg4;
   for(int i = 0; i < 12; i++){
       mtdLpii = (i+1)*5;
       sprintf(nameMtdLpii, "hmtdLpiscan_%d", mtdLpii);
       hmtdLpiscan[i] = new TH1F("hmtdLpiscan", nameMtdLpii, 100, 1330, 1480);
   }                      

   //read data
   TChain tree("T");
   //exp
   for(int i = 106; i < 128; i++){
       char inName[100];
       printf("apr07_day_%d.root\n",i);
       sprintf(inName, "/u/jkubos/analiza/gitdir/Sigma_JK/out_packed/apr07_day_%d.root", i);
//     sprintf(inName, "/home/joanna/Dokumenty/HADES/Sigma/out_packed/apr07_day_%d.root", i);
       //     sprintf(inName, "/mnt/disk1/hades/analiza/pp35/Sigma1385/pp35/exp/out_packed/apr07_day_%d.root", i);
       //     sprintf(inName, "/u/jkubos/analiza/gitdir/Sigma_JK/_exp__.root");

   //read no packed files
   //sim
   /* ifstream f("inList");
   if(!f) cout << "FAILED TO OPEN DST FILE!" << endl;
   const string line;
   while(std::getline(f, line)){
       char inName[line.size() + 1];
       strcpy(inName, line.c_str());
     cout << "_" ;
*/   tree.Add(inName);
       
     TFile *fin = TFile::Open(inName, "READ");

     TH2I *hglobal_dEdx_Mdc = (TH2I*)fin -> Get("h_ds_global_dEdx_MDC");
     TH2I *hglobal_dEdx_Mdc_acc = (TH2I*)fin -> Get("h_ds_acc_global_dEdx_MDC");

     hdEdx_Mdc -> Add(hglobal_dEdx_Mdc);
     hdEdx_Mdc_acc -> Add(hglobal_dEdx_Mdc_acc);

     fin -> Close();
   }
   cout << endl;
   
   //read info from tree
   Long_t event;
   Long_t nentries = tree.GetEntries();
   Float_t Lm, pm, pimm, pipm, Sm, 
     Lth, pth, pimth, pipth, Sth,
     Lp, pp, pimp, pipp, Sp,
     primVertR, LDecVertR, primVertZ, LDecVertZ,
     mtd, mtdLpi,
     mmabc, mmabcpip, mmabcpim, invmac, invmab, mmab, mmbc;

   tree.SetBranchAddress("fLambda_M", &Lm);
   tree.SetBranchAddress("fpartA_M", &pm);
   tree.SetBranchAddress("fpartB_M", &pimm);
   tree.SetBranchAddress("fpartC_M", &pipm);
   tree.SetBranchAddress("fSigma_M", &Sm);
   tree.SetBranchAddress("fLambda_P", &Lp);
   tree.SetBranchAddress("fpartA_P", &pp);
   tree.SetBranchAddress("fpartB_P", &pimp);
   tree.SetBranchAddress("fpartC_P", &pipp);
   tree.SetBranchAddress("fSigma_P", &Sp);
   tree.SetBranchAddress("fPrimaryVertexR", &primVertR);
   tree.SetBranchAddress("fLambdaDecayR", &LDecVertR);
   tree.SetBranchAddress("fPrimaryVertexZ", &primVertZ);
   tree.SetBranchAddress("fLambdaDecayZ", &LDecVertZ);
   tree.SetBranchAddress("fMinTrackDist", &mtd);
   tree.SetBranchAddress("fMinTrackDistLpi", &mtdLpi);
   tree.SetBranchAddress("fMMABC", &mmabc);
   tree.SetBranchAddress("fMMABCpip", &mmabcpip);
   tree.SetBranchAddress("fMMABCpim", &mmabcpim);
   tree.SetBranchAddress("finvMAC", &invmac);
   tree.SetBranchAddress("finvMAB", &invmab);
   tree.SetBranchAddress("fMMAB", &mmab);
   tree.SetBranchAddress("fMMBC", &mmbc);

   int eventNo = -1;
   TH1I *heventNo = new TH1I();
   heventNo -> SetName("heventNo");
   //event loop
   for (event = 0; event < nentries; event++){
     if(!(event%10000)){
       printf("Event: %ld\r",event);
     }
     tree.GetEntry(event); 

     //cut on miss mass ABC
     if(mmabc < cmm1 || mmabcpim < cmm2)
	 continue;
          
     //no cuts
     hLm -> Fill(Lm);
     hLp -> Fill(Lp);
     hpm -> Fill(pm);
     hpp -> Fill(pp);
     hpimm -> Fill(pimm);
     hpimp -> Fill(pimp);
     hpipm -> Fill(pipm);                                                       
     hpipp -> Fill(pipp);   
     hSm -> Fill(Sm);
     hSp -> Fill(Sp);
     hPrimVertZ -> Fill(primVertZ);
     hLDecVertZ -> Fill(LDecVertZ);
     hPrimVertR -> Fill(primVertR);
     hLDecVertR -> Fill(LDecVertR);
     hdLVR -> Fill(LDecVertR - primVertR);
     hMinTrackDist -> Fill(mtd);
     hMinTrackDistLpi -> Fill(mtdLpi);
     hMMABC -> Fill(mmabc);
     hMMABCpip -> Fill(mmabcpip);
     hMMABCpim -> Fill(mmabcpim);
     hinvMAC -> Fill(invmac);
     hinvMAB -> Fill(invmab);
     hMMAB -> Fill(mmab);
     hMMBC -> Fill(mmbc);
     hDN -> Fill(invmac, mmabcpim);

     /*   //MTD L scan
     for(int i = 0; i < 15; i++){
	 mtdLi = (i+1)*2;
	 if(mtd < mtdLi){
	     hmtdLscan[i] -> Fill(Lm);
	     //
	     for(int j = 0; j < 15; j++){
		 dLVi_mtdL = (j+1)*5;
		 if((LDecVertR - primVertR) > dLVi_mtdL)// && (LDecVertR - primVertR) < csVL2)
		     hdLVmtdLscan[i][j] -> Fill(Lm);
	     }
	 }
     }	 
     //dLV scan
     for(int i = 0; i < 15; i++){
	 dLVi = (i+1)*5;
	 if((LDecVertR - primVertR) > dLVi && (LDecVertR - primVertR) < csVL2)
	     hdLVscan[i] -> Fill(Lm);
     }	 
     */    
     //MTD S scan
     if(mtd < cmtd && (LDecVertR - primVertR) > csVL1 && Lm > cLm1 && Lm < cLm2){
       for(int i = 0; i < 12; i++){
	   mtdLpii = (i+1)*5;
	   if(mtdLpi < mtdLpii)
	       hmtdLpiscan[i] -> Fill(Sm);
       }
     }
        
     //L & S inv mass histos
     //MTD L
     if(mtd < cmtd){
       hLm_mtdL -> Fill(Lm);
       hSm_mtdL -> Fill(Sm);
       //vertLz
       //if(LDecVertZ > cvertLz1 && LDecVertZ < cvertLz2){
       //	 hLm_mtdL_vertLz -> Fill(Lm);
       //	 hSm_mtdL_vertLz -> Fill(Sm);
       //LDecVertR-PrimVertR
       if((LDecVertR - primVertR) > csVL1){// && (LDecVertR - primVertR) < csVL2){
	   hLm_mtdL_dLV -> Fill(Lm);
	   hSm_mtdL_dLV -> Fill(Sm);
	   if(Lm > cLm1 && Lm < cLm2){
	     hLm_mtdL_dLV_Lm -> Fill(Lm);
	     hSm_mtdL_dLV_Lm -> Fill(Sm);
	     if(mtdLpi < cmtdLpi){
		 hSm_mtdS -> Fill(Sm);
		 if(primVertZ > cpvz1 && primVertZ < cpvz2)
		     hSm_mtdS_pvz -> Fill(Sm);
	     }
	   }
         }
	 //}
     }
     //dLV
     if((LDecVertR - primVertR) > csVL1){// && (LDecVertR - primVertR) < csVL2){
	 hLm_dLV -> Fill(Lm);
	 hSm_dLV -> Fill(Sm);
	 if(Lm > cLm1 && Lm < cLm2){
	     hLm_dLV_Lm -> Fill(Lm);
	     hSm_dLV_Lm -> Fill(Sm);
	 }
     }
     eventNo = event;
     heventNo -> Fill(eventNo);
   }//end event loop 
   printf("\n");

   //MTD L scan
/*   for(int i = 0; i < 15; i++){
       mtdLi = (i+1)*2;
       //calc Sgn for each mtdL value
       cMtdScan -> cd(i+1);
       TH1F *htmp = (TH1F*)hmtdLscan[i] -> Clone("htmp");
       const char *hnametmp = htmp -> GetName();
       char namePeaktmp[50];
       sprintf(namePeaktmp, "%s_peak", hnametmp);
       
       cout << "--" << i << "--" << endl;
       TF1 * ftmp = new TF1("ftmp", "gaus(0)+gaus(3)+pol2(6)", 1090, 1140);
       ftmp -> SetParameters(
	   3.92029e+03,1.11517e+03,1.64829e+00,
	   1.30091e+03,1.11477e+03,4.04957e+00,
	   -1.04526e+06,9.41713e+02
	   );
       ftmp -> SetParLimits(0, 0, 100000);
       ftmp -> SetParLimits(1, 1110, 1120);
       ftmp -> SetParLimits(2, 0, 10);
       ftmp -> SetParLimits(3, 0, 100000);
       ftmp -> SetParLimits(4, 1080, 1120);
       ftmp -> SetParLimits(5, 0, 10);
       htmp -> Fit("ftmp", "0", "", 1100, 1130);
       htmp -> Fit("ftmp", "0", "", 1100, 1130);
       htmp -> SetName(namePeaktmp);
       TF1 * fsigtmp = new TF1("fsigtmp", "gaus(0)+gaus(3)", 1105, 1125);
       TF1 * fbgtmp = new TF1("fbgtmp", "pol2(0)", 1105, 1125);
       double partmp[12];
       ftmp -> GetParameters(partmp);
       fsigtmp -> SetParameters(partmp);
       fbgtmp -> SetParameters(&partmp[6]);
       
       TH1F *histSigtmp = (TH1F*)htmp -> Clone("histSigtmp");
       TH1F *histBGtmp = (TH1F*)htmp -> Clone("histBGtmp");
       histSigtmp -> Add(fbgtmp, -1);
       histSigtmp -> SetLineColor(kGreen);
       histBGtmp -> Add(fsigtmp, -1);
       histBGtmp -> SetLineColor(kRed);
       histBGtmp -> SetMarkerColor(kRed);
       histBGtmp -> SetMarkerStyle(20);
       histBGtmp -> SetMarkerSize(.7);
       
       double as = histSigtmp -> FindBin(1110);
       double bs = histSigtmp -> FindBin(1120);
       double ab = histBGtmp -> FindBin(1110);
       double bb = histBGtmp -> FindBin(1120);
       sg = histSigtmp -> Integral(as,bs);
       bg = histBGtmp -> Integral(ab,bb);
       histSigtmp -> GetXaxis() -> SetRangeUser(1110,1120);
       histBGtmp -> GetXaxis() -> SetRangeUser(1090,1140);
       
       htmp -> Draw();
       histBGtmp -> Draw("same p");
       histSigtmp -> Draw("same");
       
       mtdVal[i] = mtdLi;
       sgnVal[i] = sg/sqrt((sg*sg)+(bg*bg));
   }

   //dLV scan
   for(int i = 0; i < 15; i++){
       dLVi = (i+1)*5;
       //calc Sgn for each dLV value
       cdLVScan -> cd(i+1);
       TH1F *htmp = (TH1F*)hdLVscan[i] -> Clone("htmp");
       const char *hnametmp = htmp -> GetName();
       char namePeaktmp[50];
       sprintf(namePeaktmp, "%s_peak", hnametmp);
       
       cout << "--0" << i << "--" << endl;
       if(i < 7){
	   TF1 * ftmp = new TF1("ftmp", "gaus(0)+gaus(3)+pol3(6)", 1090, 1140);
	   ftmp -> SetParameters(
	       3.92029e+03,1.11517e+03,1.64829e+00,
	       1.30091e+03,1.11477e+03,4.04957e+00,
	       -1.04526e+06,9.41713e+02
	       );
	   TF1 * fbgtmp = new TF1("fbgtmp", "pol3(0)", 1090, 1140);
       }else{
	   TF1 * ftmp = new TF1("ftmp", "gaus(0)+gaus(3)+pol1(6)", 1090, 1140);
	   ftmp -> SetParameters(
	       3.92029e+03,1.11517e+03,1.64829e+00,
	       1.30091e+03,1.11477e+03,4.04957e+00,
	       -1.04526e+06
	       );
	   TF1 * fbgtmp = new TF1("fbgtmp", "pol1(0)", 1090, 1140);
       }
       TF1 * fsigtmp = new TF1("fsigtmp", "gaus(0)+gaus(3)", 1090, 1140);
       
       ftmp -> SetParLimits(0, 0, 100000);
       ftmp -> SetParLimits(1, 1110, 1120);
       ftmp -> SetParLimits(2, 0, 10);
       ftmp -> SetParLimits(3, 0, 100000);
       ftmp -> SetParLimits(4, 1080, 1120);
       ftmp -> SetParLimits(5, 0, 10);
       htmp -> Fit("ftmp", "0", "", 1090, 1140);
       htmp -> Fit("ftmp", "0", "", 1090, 1140);
       htmp -> SetName(namePeaktmp);

       double partmp[12];
       ftmp -> GetParameters(partmp);
       fsigtmp -> SetParameters(partmp);
       fbgtmp -> SetParameters(&partmp[6]);
       
       TH1F *histSigtmp = (TH1F*)htmp -> Clone("histSigtmp");
       TH1F *histBGtmp = (TH1F*)htmp -> Clone("histBGtmp");
       histSigtmp -> Add(fbgtmp, -1);
       histSigtmp -> SetLineColor(kGreen);
       histBGtmp -> Add(fsigtmp, -1);
       histBGtmp -> SetLineColor(kRed);
       histBGtmp -> SetMarkerColor(kRed);
       histBGtmp -> SetMarkerStyle(20);
       histBGtmp -> SetMarkerSize(.7);
       
       double as = histSigtmp -> FindBin(1110);
       double bs = histSigtmp -> FindBin(1120);
       double ab = histBGtmp -> FindBin(1110);
       double bb = histBGtmp -> FindBin(1120);
       sg2 = histSigtmp -> Integral(as,bs);
       bg2 = histBGtmp -> Integral(ab,bb);
       histSigtmp -> GetXaxis() -> SetRangeUser(1110,1120);
       histBGtmp -> GetXaxis() -> SetRangeUser(1090,1140);

       htmp -> Draw();
       fbgtmp -> Draw("same");
//       histBGtmp -> Draw("same p");
       histSigtmp -> Draw("same");
       
       dLVVal[i] = dLVi;
       sgnVal2[i] = sg2/sqrt((sg2*sg2)+(bg2*bg2));
   }

   //MTD L_dLV scan
   int n = 0;
   for(int i = 0; i < 15; i++){
       mtdLi_dLV = (i+1)*2;
       for(int j = 0; j < 15; j++){
	   cmtdLdLV[i] -> cd(j+1);
	   dLVi_mtdL = (j+1)*5;
	   //calc Sgn for each mtdL value
	   char namehtmp[50];
	   sprintf(namehtmp, "%htmp_mtdi%d_dLV%d", mtdLi_dLV, dLVi_mtdL);
	   TH1F *htmp = (TH1F*)hdLVmtdLscan[i][j] -> Clone(namehtmp);
	   const char *hnametmp = htmp -> GetName();
	   char namePeaktmp[50];
	   sprintf(namePeaktmp, "%s_peak", hnametmp);
       
	   cout << "--mtdL" << i << "dLV" << j << "--" << endl;
	   TF1 * ftmp = new TF1("ftmp", "gaus(0)+gaus(3)+pol3(6)", 1090, 1160);
	   ftmp -> SetParameters(
	       3.92029e+03,1.11517e+03,1.64829e+00,
	       1.30091e+03,1.11477e+03,4.04957e+00,
	       -2.38758e+04,2.19863e+01
	       );
	   TF1 * fbgtmp = new TF1("fbgtmp", "pol3(0)", 1090, 1160);
	   TF1 * fsigtmp = new TF1("fsigtmp", "gaus(0)+gaus(3)", 1090, 1140);

	   ftmp -> SetParLimits(0, 0, 100000);
	   ftmp -> SetParLimits(1, 1110, 1120);
	   ftmp -> SetParLimits(2, 0, 10);
	   ftmp -> SetParLimits(3, 0, 100000);
	   ftmp -> SetParLimits(4, 1080, 1120);
	   ftmp -> SetParLimits(5, 0, 10);
	   htmp -> Fit("ftmp", "0", "", 1105, 1150);
	   htmp -> Fit("ftmp", "0", "", 1105, 1150);
	   htmp -> SetName(namePeaktmp);

	   double partmp[12];
	   ftmp -> GetParameters(partmp);
	   fsigtmp -> SetParameters(partmp);
	   fbgtmp -> SetParameters(&partmp[6]);
       
	   TH1F *histSigtmp = (TH1F*)htmp -> Clone("histSigtmp");
	   TH1F *histBGtmp = (TH1F*)htmp -> Clone("histBGtmp");
	   histSigtmp -> Add(fbgtmp, -1);
	   histSigtmp -> SetLineColor(kGreen);
	   histBGtmp -> Add(fsigtmp, -1);
	   histBGtmp -> SetLineColor(kRed);
	   histBGtmp -> SetMarkerColor(kRed);
	   histBGtmp -> SetMarkerStyle(20);
	   histBGtmp -> SetMarkerSize(.7);
	   
	   double as = histSigtmp -> FindBin(1110);
	   double bs = histSigtmp -> FindBin(1120);
	   double ab = histBGtmp -> FindBin(1110);
	   double bb = histBGtmp -> FindBin(1120);
	   sg3 = histSigtmp -> Integral(as,bs);
	   bg3 = histBGtmp -> Integral(ab,bb);
	   histSigtmp -> GetXaxis() -> SetRangeUser(1110,1120);
	   histBGtmp -> GetXaxis() -> SetRangeUser(1090,1140);

	   htmp -> Draw();
	   histBGtmp -> Draw("same p");
	   histSigtmp -> Draw("same");

	   if(sg3 < 0 || bg3 < 0)
	       sgnVal3[i][j] = 0;
	   else
	       sgnVal3[i][j] = sg3/sqrt((sg3*sg3)+(bg3*bg3));
	   cout << "sg3, bg3, n, mtdi, dlvj, sgn: " << sg3 << " " << bg3 << " " << n << " "<< dLVi_mtdL << " " << mtdLi_dLV << " " << sgnVal3[i][j] << endl;
	   cout << ">>>>>>>>>>" << sgnVal3[i][j] << "<<<<<<<<<<<<<<" << endl;
	   gdLVmtdL -> SetPoint(n, dLVi_mtdL, mtdLi_dLV, sgnVal3[i][j]);
	   n++;
       }
   }
*/
   //MTD S scan
   for(int i = 0; i < 12; i++){
//       hmtdLpiscan[i] -> Rebin();
       mtdLpii = (i+1)*5;
       //calc Sgn for each mtdLpi value
       cMtdLpiScan -> cd(i+1);
       TH1F *htmp = (TH1F*)hmtdLpiscan[i] -> Clone("htmp");
       const char *hnametmp = htmp -> GetName();
       char namePeaktmp[50];
       sprintf(namePeaktmp, "%s_peak", hnametmp);
       
       cout << "--" << i << "--" << endl;
       TF1 * ftmp = new TF1("ftmp", "gaus(0)+gaus(3)+pol2(6)", 1370, 1405);
       ftmp -> SetParameters(
	   3.92029e+03,1385,5,
	   1.30091e+03,1385,5,
	   8
	   );
       ftmp -> SetParLimits(0, 0, 100000);
       ftmp -> SetParLimits(1, 1378, 1398);
       ftmp -> SetParLimits(2, 0, 10);
       ftmp -> SetParLimits(3, 0, 100000);
       ftmp -> SetParLimits(4, 1370, 1405);
       ftmp -> SetParLimits(5, 0, 10);
       htmp -> Fit("ftmp", "0", "", 1370, 1405);
       htmp -> Fit("ftmp", "0", "", 1370, 1405);
       htmp -> SetName(namePeaktmp);
       TF1 * fsigtmp = new TF1("fsigtmp", "gaus(0)+gaus(3)", 1370, 1405);
       TF1 * fbgtmp = new TF1("fbgtmp", "pol2(0)", 1370, 1405);
       double partmp[12];
       ftmp -> GetParameters(partmp);
       fsigtmp -> SetParameters(partmp);
       fbgtmp -> SetParameters(&partmp[6]);
       
       TH1F *histSigtmp = (TH1F*)htmp -> Clone("histSigtmp");
       TH1F *histBGtmp = (TH1F*)htmp -> Clone("histBGtmp");
       histSigtmp -> Add(fbgtmp, -1);
       histSigtmp -> SetLineColor(kGreen);
       histBGtmp -> Add(fsigtmp, -1);
       histBGtmp -> SetLineColor(kRed);
       histBGtmp -> SetMarkerColor(kRed);
       histBGtmp -> SetMarkerStyle(20);
       histBGtmp -> SetMarkerSize(.7);
       
       double as = histSigtmp -> FindBin(1378);
       double bs = histSigtmp -> FindBin(1398);
       double ab = histBGtmp -> FindBin(1378);
       double bb = histBGtmp -> FindBin(1398);
       sg4 = histSigtmp -> Integral(as,bs);
       bg4 = histBGtmp -> Integral(ab,bb);
       histSigtmp -> GetXaxis() -> SetRangeUser(1370,1405);
       histBGtmp -> GetXaxis() -> SetRangeUser(1370,1405);

       htmp -> Draw();
       histBGtmp -> Draw("same p");
       histSigtmp -> Draw("same");
       
       mtdLpiVal[i] = mtdLpii;
       sgnVal4[i] = sg4/sqrt((sg4*sg4)+(bg4*bg4));
   }

   //
      
   /* TGraph *gmtdscan = new TGraph(15, mtdVal, sgnVal);
   gmtdscan -> GetXaxis() -> SetTitle("MTD_{p#pi^{-}} [mm]");
   gmtdscan -> GetYaxis() -> SetTitle("#alpha = S/sqrt{S+B}");
   TGraph *gdLVscan = new TGraph(15, dLVVal, sgnVal2);
   gdLVscan -> GetXaxis() -> SetTitle("LVert-PrimVert dist [mm]");
   gdLVscan -> GetYaxis() -> SetTitle("#alpha = S/sqrt{S+B}");

   gdLVmtdL -> GetXaxis() -> SetTitle("LVert-PrimVert dist [mm]");
   gdLVmtdL -> GetYaxis() -> SetTitle("MTD_{p#pi^{-}} [mm]");
   gdLVmtdL -> GetZaxis() -> SetTitle("#alpha = S/sqrt{S+B}");
   */
   TGraph *gmtdLpiscan = new TGraph(10, mtdLpiVal, sgnVal4);
   gmtdLpiscan -> GetXaxis() -> SetTitle("MTD_{#Lambda#pi^{+}} [mm]");
   gmtdLpiscan -> GetYaxis() -> SetTitle("#alpha = S/sqrt{S+B}");
   
   
   //BG subtruction
   //L no cuts
   TH1F *histLnc = (TH1F*)hLm -> Clone("histLnc");
   const char *hnameLnc = histLnc -> GetName();
   char namePeakLnc[50];
   sprintf(namePeakLnc, "%s_peak", hnameLnc);

   cout << ">>>>>>>>fit BG no cuts<<<<<<<<<" << endl;
   TF1 * fhistLnc = new TF1("fhistLnc", "gaus(0)+gaus(3)+pol5(6)", 1090, 1140);
   fhistLnc -> SetParameters(
       2.67356e+03,1.11092e+03,1.00000e+01,
       9.66980e+03,1.11501e+03,1.93216e+00,
       -2.94616e+06,1.94250e+03,1.57966e+00,7.44350e-04,-1.40722e-06
       );
   fhistLnc -> SetParLimits(0, 0, 100000);
   fhistLnc -> SetParLimits(1, 1110, 1120);
   fhistLnc -> SetParLimits(2, 0, 10);
   fhistLnc -> SetParLimits(3, 0, 100000);
   fhistLnc -> SetParLimits(4, 1080, 1120);
   fhistLnc -> SetParLimits(5, 0, 10);
   histLnc -> Fit("fhistLnc", "0", "", 1090, 1140);
   histLnc -> Fit("fhistLnc", "0", "", 1090, 1140);
   histLnc -> SetName(namePeakLnc);
   TF1 * fsigLnc = new TF1("fsigLnc", "gaus(0)+gaus(3)", 1090, 1130);
   TF1 * fbgLnc = new TF1("fbgLnc", "pol5(0)", 1090, 1130);
   double parLnc[12];
   fhistLnc -> GetParameters(parLnc);
   fsigLnc -> SetParameters(parLnc);
   fbgLnc -> SetParameters(&parLnc[6]);

   TH1F *histSigLnc = (TH1F*)histLnc -> Clone("histSigLnc");
   TH1F *histBGLnc = (TH1F*)histLnc -> Clone("histBGLnc");
   histSigLnc -> Add(fbgLnc, -1);
   histSigLnc -> SetLineColor(kGreen);
   histBGLnc -> Add(fsigLnc, -1);
   histBGLnc -> SetLineColor(kRed);
   histBGLnc -> SetMarkerColor(kRed);
   histBGLnc -> SetMarkerStyle(20);
   histBGLnc -> SetMarkerSize(.7);

   double a = histSigLnc -> FindBin(1110);
   double b = histSigLnc -> FindBin(1120);
   double cntLmSigLnc = histSigLnc -> Integral(a,b);
   double cntLmBGLnc = histBGLnc -> Integral(a,b);
   histSigLnc -> GetXaxis() -> SetRangeUser(1110,1120);
   histBGLnc -> GetXaxis() -> SetRangeUser(1090,1140);
   /////

   //L mtd
   TH1F *histLmtd = (TH1F*)hLm_mtdL -> Clone("histLmtd");
   const char *hnameLmtd = histLmtd -> GetName();
   char namePeakLmtd[50];
   sprintf(namePeakLmtd, "%s_peak", hnameLmtd);

   cout << ">>>>>>>>fit BG MtdL cut<<<<<<<<<" << endl;
   TF1 * fhistLmtd = new TF1("fhistLmtd", "gaus(0)+gaus(3)+pol5(6)", 1090, 1140);
   fhistLmtd -> SetParameters(
       7.97144e+03,1.11503e+03,1.72359e+00,
       1.71779e+03,1.11377e+03,3.67803e+00,
       -3.46965e+06,2.17130e+03,2.01357e+00,9.48373e-04,-1.76904e-06
       );
   fhistLmtd -> SetParLimits(0, 0, 100000);
   fhistLmtd -> SetParLimits(1, 1110, 1120);
   fhistLmtd -> SetParLimits(2, 0, 10);
   fhistLmtd -> SetParLimits(3, 0, 100000);
   fhistLmtd -> SetParLimits(4, 1080, 1120);
   fhistLmtd -> SetParLimits(5, 0, 10);
   histLmtd -> Fit("fhistLmtd", "0", "", 1090, 1140);
   histLmtd -> Fit("fhistLmtd", "0", "", 1090, 1140);
   histLmtd -> SetName(namePeakLmtd);
   TF1 * fsigLmtd = new TF1("fsigLmtd", "gaus(0)+gaus(3)", 1090, 1130);
   TF1 * fbgLmtd = new TF1("fbgLmtd", "pol5(0)", 1090, 1130);
   double parLmtd[12];
   fhistLmtd -> GetParameters(parLmtd);
   fsigLmtd -> SetParameters(parLmtd);
   fbgLmtd -> SetParameters(&parLmtd[6]);

   TH1F *histSigLmtd = (TH1F*)histLmtd -> Clone("histSigLmtd");
   TH1F *histBGLmtd = (TH1F*)histLmtd -> Clone("histBGLmtd");
   histSigLmtd -> Add(fbgLmtd, -1);
   histSigLmtd -> SetLineColor(kGreen);
   histBGLmtd -> Add(fsigLmtd, -1);
   histBGLmtd -> SetLineColor(kRed);
   histBGLmtd -> SetMarkerColor(kRed);
   histBGLmtd -> SetMarkerStyle(20);
   histBGLmtd -> SetMarkerSize(.7);

   double a = histSigLmtd -> FindBin(1110);
   double b = histSigLmtd -> FindBin(1120);
   double cntLmSigLmtd = histSigLmtd -> Integral(a,b);
   double cntLmBGLmtd = histBGLmtd -> Integral(a,b);
   histSigLmtd -> GetXaxis() -> SetRangeUser(1110,1120);
   histBGLmtd -> GetXaxis() -> SetRangeUser(1090,1140);
   /////

   //L mtd dLV
   TH1F *histLmtdDlv = (TH1F*)hLm_mtdL_dLV -> Clone("histLmtdDlv");
   const char *hnameLmtdDlv = histLmtdDlv -> GetName();
   char namePeakLmtdDlv[50];
   sprintf(namePeakLmtdDlv, "%s_peak", hnameLmtdDlv);

   cout << ">>>>>>>>fit BG MtdL dLV cuts<<<<<<<<<" << endl;
   TF1 * fhistLmtdDlv = new TF1("fhistLmtdDlv", "gaus(0)+gaus(3)+pol5(6)", 1090, 1140);
   fhistLmtdDlv -> SetParameters(
			  3.92029e+03,1.11517e+03,1.64829e+00,
			  1.30091e+03,1.11477e+03,4.04957e+00,
			  -1.04526e+06,9.41713e+02,7.93433e-01,-7.10776e-04
   			  );
   fhistLmtdDlv -> SetParLimits(0, 0, 100000);
   fhistLmtdDlv -> SetParLimits(1, 1110, 1120);
   fhistLmtdDlv -> SetParLimits(2, 0, 10);
   fhistLmtdDlv -> SetParLimits(3, 0, 100000);
   fhistLmtdDlv -> SetParLimits(4, 1080, 1120);
   fhistLmtdDlv -> SetParLimits(5, 0, 10);
   histLmtdDlv -> Fit("fhistLmtdDlv", "0", "", 1090, 1140);
   histLmtdDlv -> Fit("fhistLmtdDlv", "0", "", 1090, 1140);
   histLmtdDlv -> SetName(namePeakLmtdDlv);
   TF1 * fsigLmtdDlv = new TF1("fsigLmtdDlv", "gaus(0)+gaus(3)", 1105, 1125);
   TF1 * fbgLmtdDlv = new TF1("fbgLmtdDlv", "pol5(0)", 1090, 1140);
   double parLmtdDlv[12];
   fhistLmtdDlv -> GetParameters(parLmtdDlv);
   fsigLmtdDlv -> SetParameters(parLmtdDlv);
   fbgLmtdDlv -> SetParameters(&parLmtdDlv[6]);

   TH1F *histSigLmtdDlv = (TH1F*)histLmtdDlv -> Clone("histSigLmtdDlv");
   TH1F *histBGLmtdDlv = (TH1F*)histLmtdDlv -> Clone("histBGLmtdDlv");
   histSigLmtdDlv -> Add(fbgLmtdDlv, -1);
   histSigLmtdDlv -> SetLineColor(kGreen);
   histBGLmtdDlv -> Add(fsigLmtdDlv, -1);
   histBGLmtdDlv -> SetLineColor(kRed);
   histBGLmtdDlv -> SetMarkerColor(kRed);
   histBGLmtdDlv -> SetMarkerStyle(20);
   histBGLmtdDlv -> SetMarkerSize(.7);

   double a = histSigLmtdDlv -> FindBin(1110);
   double b = histSigLmtdDlv -> FindBin(1120);
   double cntLmSigLmtdDlv = histSigLmtdDlv -> Integral(a,b);
   double cntLmBGLmtdDlv = histBGLmtdDlv -> Integral(a,b);
   histSigLmtdDlv -> GetXaxis() -> SetRangeUser(1110,1120);
   histBGLmtdDlv -> GetXaxis() -> SetRangeUser(1090,1140);
   ///
   //L dLV cut
   TH1F *histLDlv = (TH1F*)hLm_dLV -> Clone("histLDlv");
   const char *hnameLDlv = histLDlv -> GetName();
   char namePeakLDlv[50];
   sprintf(namePeakLDlv, "%s_peak", hnameLDlv);

   cout << ">>>>>>>>fit BG dLV cuts<<<<<<<<<" << endl;
   TF1 * fhistLDlv = new TF1("fhistLDlv", "gaus(0)+gaus(3)+pol6(6)", 1090, 1140);
   fhistLDlv -> SetParameters(
       3.22861e+03,1.11505e+03,2.26541e+00,
       4.49200e+02,1.10917e+03,1.79182e+00,
       -5.84585e+05,8.67222e+02,5.07379e-01,-9.69784e-04,-6.22429e-08
       );
   fhistLDlv -> SetParLimits(0, 0, 100000);
   fhistLDlv -> SetParLimits(1, 1110, 1120);
   fhistLDlv -> SetParLimits(2, 0, 10);
   fhistLDlv -> SetParLimits(3, 0, 100000);
   fhistLDlv -> SetParLimits(4, 1080, 1120);
   fhistLDlv -> SetParLimits(5, 0, 10);
   histLDlv -> Fit("fhistLDlv", "0", "", 1090, 1140);
   histLDlv -> Fit("fhistLDlv", "0", "", 1090, 1140);
   histLDlv -> SetName(namePeakLDlv);
   TF1 * fsigLDlv = new TF1("fsigLDlv", "gaus(0)+gaus(3)", 1105, 1125);
   TF1 * fbgLDlv = new TF1("fbgLDlv", "pol6(0)", 1090, 1140);
   double parLDlv[12];
   fhistLDlv -> GetParameters(parLDlv);
   fsigLDlv -> SetParameters(parLDlv);
   fbgLDlv -> SetParameters(&parLDlv[6]);

   TH1F *histSigLDlv = (TH1F*)histLDlv -> Clone("histSigLDlv");
   TH1F *histBGLDlv = (TH1F*)histLDlv -> Clone("histBGLDlv");
   histSigLDlv -> Add(fbgLDlv, -1);
   histSigLDlv -> SetLineColor(kGreen);
   histBGLDlv -> Add(fsigLDlv, -1);
   histBGLDlv -> SetLineColor(kRed);
   histBGLDlv -> SetMarkerColor(kRed);
   histBGLDlv -> SetMarkerStyle(20);
   histBGLDlv -> SetMarkerSize(.7);

   double a = histSigLDlv -> FindBin(1110);
   double b = histSigLDlv -> FindBin(1120);
   double cntLmSigLDlv = histSigLDlv -> Integral(a,b);
   double cntLmBGLDlv = histBGLDlv -> Integral(a,b);
   histSigLDlv -> GetXaxis() -> SetRangeUser(1110,1120);
   histBGLDlv -> GetXaxis() -> SetRangeUser(1090,1140);
   ///
   //>>>>>>>>>>>>>Sigma<<<<<<<<<<<<<<<<<<<<<<<
   /*  //S no cuts
   TH1F *histSnc = (TH1F*)hSm -> Clone("histSnc");
   const char *hnameSnc = histSnc -> GetName();
   char namePeakSnc[50];
   sprintf(namePeakSnc, "%s_peak", hnameSnc);

   cout << ">>>>>>>>fit BG no cuts<<<<<<<<<" << endl;
   TF1 * fhistSnc = new TF1("fhistSnc", "gaus(0)+gaus(3)+pol4(6)", 1090, 1140);
   fhistSnc -> SetParameters(
       3.92029e+03,1.11517e+03,1.64829e+00,
       1.30091e+03,1.11477e+03,4.04957e+00,
       -1.04526e+06,9.41713e+02,7.93433e-01,-7.10776e-04
       );
   fhistSnc -> SetParLimits(0, 0, 100000);
   fhistSnc -> SetParLimits(1, 1110, 1120);
   fhistSnc -> SetParLimits(2, 0, 10);
   fhistSnc -> SetParLimits(3, 0, 100000);
   fhistSnc -> SetParLimits(4, 1080, 1120);
   fhistSnc -> SetParLimits(5, 0, 10);
   histSnc -> Fit("fhistSnc", "0", "", 1090, 1140);
   histSnc -> Fit("fhistSnc", "0", "", 1090, 1140);
   histSnc -> SetName(namePeakSnc);
   TF1 * fsigSnc = new TF1("fsigSnc", "gaus(0)+gaus(3)", 1105, 1125);
   TF1 * fbgSnc = new TF1("fbgSnc", "pol4(0)", 1090, 1140);
   double parSnc[12];
   fhistSnc -> GetParameters(parSnc);
   fsigSnc -> SetParameters(parSnc);
   fbgSnc -> SetParameters(&parSnc[6]);

   TH1F *histSigSnc = (TH1F*)histSnc -> Clone("histSigSnc");
   TH1F *histBGSnc = (TH1F*)histSnc -> Clone("histBGSnc");
   histSigSnc -> Add(fbgSnc, -1);
   histSigSnc -> SetLineColor(kGreen);
   histBGSnc -> Add(fsigSnc, -1);
   histBGSnc -> SetLineColor(kRed);
   histBGSnc -> SetMarkerColor(kRed);
   histBGSnc -> SetMarkerStyle(20);
   histBGSnc -> SetMarkerSize(.7);

   double a = histSigSnc -> FindBin(1110);
   double b = histSigSnc -> FindBin(1120);
   double cntSmSigSnc = histSigSnc -> Integral(a,b);
   double cntSmBGSnc = histBGSnc -> Integral(a,b);
   histSigSnc -> GetXaxis() -> SetRangeUser(1110,1120);
   histBGSnc -> GetXaxis() -> SetRangeUser(1090,1140);
   /////
   */
   //S mtdL
/*   TH1F *histSmtdL = (TH1F*)hSm_mtdL -> Clone("histSmtdL");
   const char *hnameSmtdL = histSmtdL -> GetName();
   char namePeakSmtdL[50];
   sprintf(namePeakSmtdL, "%s_peak", hnameSmtdL);

   cout << ">>>>>>>>fit BG MtdL cut<<<<<<<<<" << endl;
   TF1 * fhistSmtdL = new TF1("fhistSmtdL", "gaus(0)+gaus(3)+pol2(6)", 1370, 1405);
   fhistSmtdL -> SetParameters(
       3.92029e+03,1385,5,
       1.30091e+03,1385,5,
       8
       );
   fhistSmtdL -> SetParLimits(0, 0, 100000);
   fhistSmtdL -> SetParLimits(1, 1378, 1398);
   fhistSmtdL -> SetParLimits(2, 0, 10);
   fhistSmtdL -> SetParLimits(3, 0, 100000);
   fhistSmtdL -> SetParLimits(4, 1370, 1405);
   fhistSmtdL -> SetParLimits(5, 0, 10);
   histSmtdL -> Fit("fhistSmtdL", "0", "", 1370, 1405);
   histSmtdL -> Fit("fhistSmtdL", "0", "", 1370, 1405);
   histSmtdL -> SetName(namePeakSmtdL);
   TF1 * fsigSmtdL = new TF1("fsigSmtdL", "gaus(0)+gaus(3)", 1370, 1405);
   TF1 * fbgSmtdL = new TF1("fbgSmtdL", "pol4(0)", 1370, 1405);
   double parSmtdL[12];
   fhistSmtdL -> GetParameters(parSmtdL);
   fsigSmtdL -> SetParameters(parSmtdL);
   fbgSmtdL -> SetParameters(&parSmtdL[6]);

   TH1F *histSigSmtdL = (TH1F*)histSmtdL -> Clone("histSigSmtdL");
   TH1F *histBGSmtdL = (TH1F*)histSmtdL -> Clone("histBGSmtdL");
   histSigSmtdL -> Add(fbgSmtdL, -1);
   histSigSmtdL -> SetLineColor(kGreen);
   histBGSmtdL -> Add(fsigSmtdL, -1);
   histBGSmtdL -> SetLineColor(kRed);
   histBGSmtdL -> SetMarkerColor(kRed);
   histBGSmtdL -> SetMarkerStyle(20);
   histBGSmtdL -> SetMarkerSize(.7);

   double a = histSigSmtdL -> FindBin(1378);
   double b = histSigSmtdL -> FindBin(1398);
   double cntSmSigSmtdL = histSigSmtdL -> Integral(a,b);
   double cntSmBGSmtdL = histBGSmtdL -> Integral(a,b);
   histSigSmtdL -> GetXaxis() -> SetRangeUser(1370,1405);
   histBGSmtdL -> GetXaxis() -> SetRangeUser(1370,1405);
   /////
   */ 
   /* //S mtdL dLV
   TH1F *histSmtdLDlv = (TH1F*)hSm_mtdL_dLV -> Clone("histSmtdLDlv");
   const char *hnameSmtdLDlv = histSmtdLDlv -> GetName();
   char namePeakSmtdLDlv[50];
   sprintf(namePeakSmtdLDlv, "%s_peak", hnameSmtdLDlv);

   cout << ">>>>>>>>fit BG MtdL dLV cuts<<<<<<<<<" << endl;
   TF1 * fhistSmtdLDlv = new TF1("fhistSmtdLDlv", "gaus(0)+gaus(3)+pol4(6)", 1090, 1140);
   fhistSmtdLDlv -> SetParameters(
			  3.92029e+03,1.11517e+03,1.64829e+00,
			  1.30091e+03,1.11477e+03,4.04957e+00,
			  -1.04526e+06,9.41713e+02,7.93433e-01,-7.10776e-04
   			  );
   fhistSmtdLDlv -> SetParLimits(0, 0, 100000);
   fhistSmtdLDlv -> SetParLimits(1, 1110, 1120);
   fhistSmtdLDlv -> SetParLimits(2, 0, 10);
   fhistSmtdLDlv -> SetParLimits(3, 0, 100000);
   fhistSmtdLDlv -> SetParLimits(4, 1080, 1120);
   fhistSmtdLDlv -> SetParLimits(5, 0, 10);
   histSmtdLDlv -> Fit("fhistSmtdLDlv", "0", "", 1090, 1140);
   histSmtdLDlv -> Fit("fhistSmtdLDlv", "0", "", 1090, 1140);
   histSmtdLDlv -> SetName(namePeakSmtdLDlv);
   TF1 * fsigSmtdLDlv = new TF1("fsigSmtdLDlv", "gaus(0)+gaus(3)", 1105, 1125);
   TF1 * fbgSmtdLDlv = new TF1("fbgSmtdLDlv", "pol4(0)", 1090, 1140);
   double parSmtdLDlv[12];
   fhistSmtdLDlv -> GetParameters(parSmtdLDlv);
   fsigSmtdLDlv -> SetParameters(parSmtdLDlv);
   fbgSmtdLDlv -> SetParameters(&parSmtdLDlv[6]);

   TH1F *histSigSmtdLDlv = (TH1F*)histSmtdLDlv -> Clone("histSigSmtdLDlv");
   TH1F *histBGSmtdLDlv = (TH1F*)histSmtdLDlv -> Clone("histBGSmtdLDlv");
   histSigSmtdLDlv -> Add(fbgSmtdLDlv, -1);
   histSigSmtdLDlv -> SetLineColor(kGreen);
   histBGSmtdLDlv -> Add(fsigSmtdLDlv, -1);
   histBGSmtdLDlv -> SetLineColor(kRed);
   histBGSmtdLDlv -> SetMarkerColor(kRed);
   histBGSmtdLDlv -> SetMarkerStyle(20);
   histBGSmtdLDlv -> SetMarkerSize(.7);

   double a = histSigSmtdLDlv -> FindBin(1110);
   double b = histSigSmtdLDlv -> FindBin(1120);
   double cntSmSigSmtdLDlv = histSigSmtdLDlv -> Integral(a,b);
   double cntSmBGSmtdLDlv = histBGSmtdLDlv -> Integral(a,b);
   histSigSmtdLDlv -> GetXaxis() -> SetRangeUser(1110,1120);
   histBGSmtdLDlv -> GetXaxis() -> SetRangeUser(1090,1140);
   ///
   *//* //S dLV cut
   TH1F *histSDlv = (TH1F*)hSm_dLV -> Clone("histSDlv");
   const char *hnameSDlv = histSDlv -> GetName();
   char namePeakSDlv[50];
   sprintf(namePeakSDlv, "%s_peak", hnameSDlv);

   cout << ">>>>>>>>fit BG dLV cuts<<<<<<<<<" << endl;
   TF1 * fhistSDlv = new TF1("fhistSDlv", "gaus(0)+gaus(3)+pol4(6)", 1090, 1140);
   fhistSDlv -> SetParameters(
       3.92029e+03,1.11517e+03,1.64829e+00,
       1.30091e+03,1.11477e+03,4.04957e+00,
       -1.04526e+06,9.41713e+02,7.93433e-01,-7.10776e-04
       );
   fhistSDlv -> SetParLimits(0, 0, 100000);
   fhistSDlv -> SetParLimits(1, 1110, 1120);
   fhistSDlv -> SetParLimits(2, 0, 10);
   fhistSDlv -> SetParLimits(3, 0, 100000);
   fhistSDlv -> SetParLimits(4, 1080, 1120);
   fhistSDlv -> SetParLimits(5, 0, 10);
   histSDlv -> Fit("fhistSDlv", "0", "", 1090, 1140);
   histSDlv -> Fit("fhistSDlv", "0", "", 1090, 1140);
   histSDlv -> SetName(namePeakSDlv);
   TF1 * fsigSDlv = new TF1("fsigSDlv", "gaus(0)+gaus(3)", 1105, 1125);
   TF1 * fbgSDlv = new TF1("fbgSDlv", "pol4(0)", 1090, 1140);
   double parSDlv[12];
   fhistSDlv -> GetParameters(parSDlv);
   fsigSDlv -> SetParameters(parSDlv);
   fbgSDlv -> SetParameters(&parSDlv[6]);

   TH1F *histSigSDlv = (TH1F*)histSDlv -> Clone("histSigSDlv");
   TH1F *histBGSDlv = (TH1F*)histSDlv -> Clone("histBGSDlv");
   histSigSDlv -> Add(fbgSDlv, -1);
   histSigSDlv -> SetLineColor(kGreen);
   histBGSDlv -> Add(fsigSDlv, -1);
   histBGSDlv -> SetLineColor(kRed);
   histBGSDlv -> SetMarkerColor(kRed);
   histBGSDlv -> SetMarkerStyle(20);
   histBGSDlv -> SetMarkerSize(.7);

   double a = histSigSDlv -> FindBin(1110);
   double b = histSigSDlv -> FindBin(1120);
   double cntSmSigSDlv = histSigSDlv -> Integral(a,b);
   double cntSmBGSDlv = histBGSDlv -> Integral(a,b);
   histSigSDlv -> GetXaxis() -> SetRangeUser(1110,1120);
   histBGSDlv -> GetXaxis() -> SetRangeUser(1090,1140);
     */ ///
   //S mtdL dLV Lm
   TH1F *histSmtdLDlvLm = (TH1F*)hSm_mtdL_dLV_Lm -> Clone("histSmtdLDlvLm");
   const char *hnameSmtdLDlvLm = histSmtdLDlvLm -> GetName();
   char namePeakSmtdLDlvLm[50];
   sprintf(namePeakSmtdLDlvLm, "%s_peak", hnameSmtdLDlvLm);

   cout << ">>>>>>>>fit BG MtdL dLV Lm cuts<<<<<<<<<" << endl;
   TF1 * fhistSmtdLDlvLm = new TF1("fhistSmtdLDlvLm", "gaus(0)+gaus(3)+pol3(6)", 1340, 1395);
   fhistSmtdLDlvLm -> SetParameters(
       9.97029e+04,1.37800e+03,1.84764e-01,
       8.37420e+01,1.38182e+03,5.71984e+00,
       1.04313e+02,-5.10370e-01,5.41819e-04
       );
   fhistSmtdLDlvLm -> SetParLimits(0, 0, 100000);
   fhistSmtdLDlvLm -> SetParLimits(1, 1378, 1392);
   fhistSmtdLDlvLm -> SetParLimits(2, 0, 10);
   fhistSmtdLDlvLm -> SetParLimits(3, 0, 100000);
   fhistSmtdLDlvLm -> SetParLimits(4, 1378, 1392);
   fhistSmtdLDlvLm -> SetParLimits(5, 0, 10);
   histSmtdLDlvLm -> Fit("fhistSmtdLDlvLm", "0", "", 1370, 1395);
   histSmtdLDlvLm -> Fit("fhistSmtdLDlvLm", "0", "", 1370, 1395);
   histSmtdLDlvLm -> SetName(namePeakSmtdLDlvLm);
   TF1 * fsigSmtdLDlvLm = new TF1("fsigSmtdLDlvLm", "gaus(0)+gaus(3)", 1367, 1405);
   TF1 * fbgSmtdLDlvLm = new TF1("fbgSmtdLDlvLm", "pol3(0)", 1335, 1405);
   double parSmtdLDlvLm[12];
   fhistSmtdLDlvLm -> GetParameters(parSmtdLDlvLm);
   fsigSmtdLDlvLm -> SetParameters(parSmtdLDlvLm);
   fbgSmtdLDlvLm -> SetParameters(&parSmtdLDlvLm[6]);

   TH1F *histSigSmtdLDlvLm = (TH1F*)histSmtdLDlvLm -> Clone("histSigSmtdLDlvLm");
   TH1F *histBGSmtdLDlvLm = (TH1F*)histSmtdLDlvLm -> Clone("histBGSmtdLDlvLm");
   histSigSmtdLDlvLm -> Add(fbgSmtdLDlvLm, -1);
   histSigSmtdLDlvLm -> SetLineColor(kGreen);
   histBGSmtdLDlvLm -> Add(fsigSmtdLDlvLm, -1);
   histBGSmtdLDlvLm -> SetLineColor(kRed);
   histBGSmtdLDlvLm -> SetMarkerColor(kRed);
   histBGSmtdLDlvLm -> SetMarkerStyle(20);
   histBGSmtdLDlvLm -> SetMarkerSize(.7);

   double a = histSigSmtdLDlvLm -> FindBin(1378);
   double b = histSigSmtdLDlvLm -> FindBin(1398);
   double cntSmSigSmtdLDlvLm = histSigSmtdLDlvLm -> Integral(a,b);
   double cntSmBGSmtdLDlvLm = histBGSmtdLDlvLm -> Integral(a,b);
   histSigSmtdLDlvLm -> GetXaxis() -> SetRangeUser(1370,1405);
   histBGSmtdLDlvLm -> GetXaxis() -> SetRangeUser(1370,1405);
   ///
   //S mtdL dLV Lm mtdS
   TH1F *histSmtdS = (TH1F*)hSm_mtdS -> Clone("histSmtdS");
   const char *hnameSmtdS = histSmtdS -> GetName();
   char namePeakSmtdS[50];
   sprintf(namePeakSmtdS, "%s_peak", hnameSmtdS);

   cout << ">>>>>>>>fit BG MtdL dLV cuts<<<<<<<<<" << endl;
   TF1 * fhistSmtdS = new TF1("fhistSmtdS", "gaus(0)+gaus(3)+pol4(6)", 1330, 1405);
   fhistSmtdS -> SetParameters(
			  3.92029e+03,1.11517e+03,1.64829e+00,
			  1.30091e+03,1.11477e+03,4.04957e+00,
			  -1.04526e+06,9.41713e+02,7.93433e-01,-7.10776e-04
   			  );
   fhistSmtdS -> SetParLimits(0, 0, 100000);
   fhistSmtdS -> SetParLimits(1, 1378, 1398);
   fhistSmtdS -> SetParLimits(2, 0, 10);
   fhistSmtdS -> SetParLimits(3, 0, 100000);
   fhistSmtdS -> SetParLimits(4, 1370, 1400);
   fhistSmtdS -> SetParLimits(5, 0, 10);
   histSmtdS -> Fit("fhistSmtdS", "0", "", 1370, 1395);
   histSmtdS -> Fit("fhistSmtdS", "0", "", 1370, 1395);
   histSmtdS -> SetName(namePeakSmtdS);
   TF1 * fsigSmtdS = new TF1("fsigSmtdS", "gaus(0)+gaus(3)", 1370, 1395);
   TF1 * fbgSmtdS = new TF1("fbgSmtdS", "pol4(0)", 1360, 1405);
   double parSmtdS[12];
   fhistSmtdS -> GetParameters(parSmtdS);
   fsigSmtdS -> SetParameters(parSmtdS);
   fbgSmtdS -> SetParameters(&parSmtdS[6]);

   TH1F *histSigSmtdS = (TH1F*)histSmtdS -> Clone("histSigSmtdS");
   TH1F *histBGSmtdS = (TH1F*)histSmtdS -> Clone("histBGSmtdS");
   histSigSmtdS -> Add(fbgSmtdS, -1);
   histSigSmtdS -> SetLineColor(kGreen);
   histBGSmtdS -> Add(fsigSmtdS, -1);
   histBGSmtdS -> SetLineColor(kRed);
   histBGSmtdS -> SetMarkerColor(kRed);
   histBGSmtdS -> SetMarkerStyle(20);
   histBGSmtdS -> SetMarkerSize(.7);

   double a = histSigSmtdS -> FindBin(1378);
   double b = histSigSmtdS -> FindBin(1398);
   double cntSmSigSmtdS = histSigSmtdS -> Integral(a,b);
   double cntSmBGSmtdS = histBGSmtdS -> Integral(a,b);
   histSigSmtdS -> GetXaxis() -> SetRangeUser(1370,1405);
   histBGSmtdS -> GetXaxis() -> SetRangeUser(1370,1405);
   ///
   //S mtdL dLV Lm mtdS pvz
   TH1F *histSmtdSpvz = (TH1F*)hSm_mtdS_pvz -> Clone("histSmtdSpvz");
   const char *hnameSmtdSpvz = histSmtdSpvz -> GetName();
   char namePeakSmtdSpvz[50];
   sprintf(namePeakSmtdSpvz, "%s_peak", hnameSmtdSpvz);

   cout << ">>>>>>>>fit BG MtdL dLV Lm MtdS DLv cuts<<<<<<<<<" << endl;
   TF1 * fhistSmtdSpvz = new TF1("fhistSmtdSpvz", "gaus(0)+gaus(3)+pol4(6)", 1330, 1405);
   fhistSmtdSpvz -> SetParameters(
       3.92029e+03,1.11517e+03,1.64829e+00,
       1.30091e+03,1.11477e+03,4.04957e+00,
       -1.04526e+06,9.41713e+02,7.93433e-01,-7.10776e-04
       );
   fhistSmtdSpvz -> SetParLimits(0, 0, 100000);
   fhistSmtdSpvz -> SetParLimits(1, 1378, 1398);
   fhistSmtdSpvz -> SetParLimits(2, 0, 10);
   fhistSmtdSpvz -> SetParLimits(3, 0, 100000);
   fhistSmtdSpvz -> SetParLimits(4, 1370, 1450);
   fhistSmtdSpvz -> SetParLimits(5, 0, 10);
   histSmtdSpvz -> Fit("fhistSmtdSpvz", "0", "", 1370, 1395);
   histSmtdSpvz -> Fit("fhistSmtdSpvz", "0", "", 1370, 1395);
   histSmtdSpvz -> SetName(namePeakSmtdSpvz);
   TF1 * fsigSmtdSpvz = new TF1("fsigSmtdSpvz", "gaus(0)+gaus(3)", 1370, 1395);
   TF1 * fbgSmtdSpvz = new TF1("fbgSmtdSpvz", "pol4(0)", 1370, 1395);
   double parSmtdSpvz[12];
   fhistSmtdSpvz -> GetParameters(parSmtdSpvz);
   fsigSmtdSpvz -> SetParameters(parSmtdSpvz);
   fbgSmtdSpvz -> SetParameters(&parSmtdSpvz[6]);

   TH1F *histSigSmtdSpvz = (TH1F*)histSmtdSpvz- > Clone("histSigSmtdSpvz");
   TH1F *histBGSmtdSpvz = (TH1F*)histSmtdSpvz -> Clone("histBGSmtdSpvz");
   histSigSmtdSpvz -> Add(fbgSmtdSpvz, -1);
   histSigSmtdSpvz -> SetLineColor(kGreen);
   histBGSmtdSpvz -> Add(fsigSmtdSpvz, -1);
   histBGSmtdSpvz -> SetLineColor(kRed);
   histBGSmtdSpvz -> SetMarkerColor(kRed);
   histBGSmtdSpvz -> SetMarkerStyle(20);
   histBGSmtdSpvz -> SetMarkerSize(.7);

   double a = histSigSmtdSpvz -> FindBin(1378);
   double b = histSigSmtdSpvz -> FindBin(1398);
   double cntSmSigSmtdSpvz = histSigSmtdSpvz -> Integral(a,b);
   double cntSmBGSmtdSpvz = histBGSmtdSpvz -> Integral(a,b);
   histSigSmtdSpvz -> GetXaxis() -> SetRangeUser(1370,1395);
   histBGSmtdSpvz -> GetXaxis() -> SetRangeUser(1370,1395);
   ///

   //end BG
   
   cout << "nentries=" << nentries << endl << "S: " << endl << "Lambda:" << " Lnc=" << cntLmSigLnc << " LmtdL=" << cntLmSigLmtd << " LmtdLdLV=" << cntLmSigLmtdDlv << " LdLV=" << cntLmSigLDlv << endl;
   //  cout << "Sigma: " << " Snc=" << cntSmSigSnc << " SmtdL=" << cntSmSigSmtdL << " SmtdLdLV=" << cntSmSigSmtdLDlv << " SdLV=" << cntSmSigSDlv << " SLcuts=" << cntSmSigSmtdLDlvLm << " SLcutsmtdS=" << cntSmSigSmtdS << endl;
   cout << "B: " << endl << "Lambda:" << " Lnc=" << cntLmBGLnc << " LmtdL=" << cntLmBGLmtd << " LmtdLdLV=" << cntLmBGLmtdDlv << " LdLV=" << cntLmBGLDlv << endl;
   //  cout << "Sigma: " << " Snc=" << cntSmBGSnc << " SmtdL=" << cntSmBGSmtdL << " SmtdLdLV=" << cntSmBGSmtdLDlv << " SdLV=" << cntSmBGSDlv << " SLcuts=" << cntSmBGSmtdLDlvLm << " SLcutsmtdS=" << cntSmBGSmtdS << endl;       


   //reco effi
   double effLS, effLBG, effSS, effSBG;
   effLS = cntLmSigLmtdDlv/nentries*100;
   effLBG = cntLmBGLmtdDlv/nentries*100;
//   effSS = cntSmSigSmtdLDlv/nentries*100;
//   effSBG = cntSmBGSmtdLDlv/nentries*100;
//   cout << "Effi: LMtdDlvLm = " << effLS << " LMtdDlvLmBG = " << effLBG << " SMtdLDlvSm = " << effSS << " SMtdLDlvSmBG = " << effSBG << endl;

     
   //edit histos
   hdEdx_Mdc -> GetXaxis() -> SetTitle("p*q [q*MeV/c]");
   hdEdx_Mdc -> GetYaxis() -> SetTitle("dEdx [a.u.]");
   hdEdx_Mdc_acc -> GetXaxis() -> SetTitle("p*q [q*MeV/c]");
   hdEdx_Mdc_acc -> GetYaxis() -> SetTitle("dEdx [a.u.]");

   hLm -> GetXaxis() -> SetTitle("m_{p#pi^{-}} [MeV]");
   hpm -> GetXaxis() -> SetTitle("m_{p} [MeV]");
   hpimm -> GetXaxis() -> SetTitle("m_{#pi^{-}} [MeV]");
   hpipm -> GetXaxis() -> SetTitle("m_{#pi^{+}} [MeV]");
   hSm -> GetXaxis() -> SetTitle("m_{p#pi^{-}#pi^{+}} [MeV]");
   hLp -> GetXaxis() -> SetTitle("p_{p#pi^{-}} [MeV/c]");
   hpp -> GetXaxis() -> SetTitle("p_{p} [MeV/c]");
   hpimp -> GetXaxis() -> SetTitle("p_{#pi^{-}} [MeV/c]");
   hpipp -> GetXaxis() -> SetTitle("p_{#pi^{+}} [MeV/c]");
   hSp -> GetXaxis() -> SetTitle("p_{p#pi^{-}#pi^{+}} [MeV/c]");

   hPrimVertZ -> GetXaxis() -> SetTitle("PrimVertZ [mm]");
   hLDecVertZ -> GetXaxis() -> SetTitle("#Lambda DecayVertZ [mm]");
   hPrimVertR -> GetXaxis() -> SetTitle("PrimVertR [mm]");
   hLDecVertR -> GetXaxis() -> SetTitle("#Lambda DecayVertR [mm]");
   hdLVR -> GetXaxis() -> SetTitle("#Lambda DecVartR - PrimVertR [mm]");
   hMinTrackDist -> GetXaxis() -> SetTitle("MTD_{p#pi^{-}} [mm]");
   hMinTrackDistLpi -> GetXaxis() -> SetTitle("MTD_{#Lambda#pi^{+}} [mm]");

   //   hLm_mtdL_vertLz -> GetXaxis() -> SetTitle("m_{p#pi^{-}} [MeV]");
   //   hSm_mtdL_vertLz -> GetXaxis() -> SetTitle("m_{p#pi^{-}#pi^{+}} [MeV]");
   hLm_mtdL -> GetXaxis() -> SetTitle("m_{p#pi^{-}} [MeV]");
   hSm_mtdL -> GetXaxis() -> SetTitle("m_{p#pi^{-}#pi^{+}} [MeV]");
   hLm_mtdL_dLV -> GetXaxis() -> SetTitle("m_{p#pi^{-}} [MeV]");
   hSm_mtdL_dLV -> GetXaxis() -> SetTitle("m_{p#pi^{-}#pi^{+}} [MeV]");
   hLm_mtdL_dLV_Lm -> GetXaxis() -> SetTitle("m_{p#pi^{-}} [MeV]");
   hSm_mtdL_dLV_Lm -> GetXaxis() -> SetTitle("m_{p#pi^{-}#pi^{+}} [MeV]");
   hLm_dLV -> GetXaxis() -> SetTitle("m_{p#pi^{-}} [MeV]");
   hSm_dLV -> GetXaxis() -> SetTitle("m_{p#pi^{-}#pi^{+}} [MeV]");
   hLm_dLV_Lm -> GetXaxis() -> SetTitle("m_{p#pi^{-}} [MeV]");
   hSm_dLV_Lm -> GetXaxis() -> SetTitle("m_{p#pi^{-}#pi^{+}} [MeV]");
   hSm_mtdS -> GetXaxis() -> SetTitle("m_{p#pi^{-}#pi^{+}} [MeV]");
   hSm_mtdS_pvz -> GetXaxis() -> SetTitle("m_{p#pi^{-}#pi^{+}} [MeV]");

   //   hLmRecBG -> GetXaxis() -> SetTitle("m_{p#pi^{-}} [MeV]");
   // histSig -> GetXaxis() -> SetTitle("m_{p#pi^{-}} [MeV]");
   // histBG -> GetXaxis() -> SetTitle("m_{p#pi^{-}} [MeV]");
   
   hMMABC -> GetXaxis() -> SetTitle("MM_{p#pi^{-}#pi^{+}}");
   hMMABCpip -> GetXaxis() -> SetTitle("MM_{p#pi^{-}#pi^{+}} + M_{#pi^{+}}");
   hMMABCpim -> GetXaxis() -> SetTitle("MM_{p#pi^{-}#pi^{+}} + M_{#pi^{-}}");
   hinvMAC -> GetXaxis() -> SetTitle("invM_{p#pi^{+}}");
   hinvMAB -> GetXaxis() -> SetTitle("invM_{p#pi^{-}}");

   hDN -> GetXaxis() -> SetTitle("invM_{p#pi^{-}}");
   hDN -> GetYaxis() -> SetTitle("MM_{p#pi^{-}#pi^{+}} + M_{#pi^{+}}");

   //graph cuts
   TCanvas *cPidCuts = new TCanvas("cPidCuts", "cPidCuts");
   cPidCuts -> cd();
   gPad -> SetLogz();
   hdEdx_Mdc -> Draw("colz");
   gCutp -> Draw("same");
   gCutpim -> Draw("same");
   gCutpip -> Draw("same");

   //BG
   //L
   cLsbg -> cd(1);
   hLm -> Draw();
   fbgLnc -> Draw("same p");
   fsigLnc -> Draw("same p");
   histSigLnc -> Draw("same");
   cLsbg -> cd(2);
   hLm_mtdL -> Draw();
   fbgLmtd -> Draw("same p");
   fsigLmtd -> Draw("same p");
   histSigLmtd -> Draw("same");
   cLsbg -> cd(3);
   hLm_dLV -> Draw();
   fbgLDlv -> Draw("same p");
   fsigLDlv -> Draw("same p");
   histSigLDlv -> Draw("same");
   cLsbg -> cd(4);
   hLm_mtdL_dLV -> Draw();
   fbgLmtdDlv -> Draw("same p");
   fsigLmtdDlv -> Draw("same p");
   histSigLmtdDlv -> Draw("same");
   //S
   cSsbg -> cd(1);
   hSm -> Draw();
   //fbgSnc -> Draw("same p");
   //histSigSnc -> Draw("same");
   cSsbg -> cd(2);
   hSm_mtdL_dLV_Lm -> Draw();
   fsigSmtdLDlvLm -> Draw("same p");
   histSigSmtdLDlvLm -> Draw("same");
   cSsbg -> cd(3);
   hSm_mtdS -> Draw();
   fsigSmtdS -> Draw("same p");
   histSigSmtdS -> Draw("same");
   cSsbg -> cd(4);
   hSm_mtdS_pvz -> Draw();
   fsigSmtdSpvz -> Draw("same p");
   histSigSmtdSpvz -> Draw("same");

   
   //writing histos
   TFile *fout = TFile::Open("./outputs/anaSigmaOut_exp_all_nmm.root", "RECREATE");
   hLm -> Write();
   hpm -> Write();
   hpimm -> Write();
   hpipm -> Write();
   hSm -> Write();
   hLp -> Write();
   hpp -> Write();
   hpimp -> Write();
   hpipp -> Write();
   hSp -> Write();
   
   hdEdx_Mdc -> Write();
   hdEdx_Mdc_acc -> Write();
   cPidCuts -> Write();
/*   cMtdScan -> Write();
   gmtdscan -> Write();
   cdLVScan -> Write();
   gdLVscan -> Write();
   gdLVmtdL -> Write();
   for(int i = 0; i < 15; i++)
       cmtdLdLV[i] -> Write();
*/   cMtdLpiScan -> Write();
   gmtdLpiscan -> Write();

   hPrimVertZ -> Write();
   hLDecVertZ -> Write();
   hPrimVertR -> Write();
   hLDecVertR -> Write();
   hdLVR -> Write();
   hMinTrackDist -> Write();
   hMinTrackDistLpi -> Write();

   hLm_mtdL -> Write();
   hSm_mtdL -> Write();
   //   hLm_mtdL_vertLz -> Write();
   //   hSm_mtdL_vertLz -> Write();
   hLm_mtdL_dLV -> Write();
   hSm_mtdL_dLV -> Write();
   hLm_mtdL_dLV_Lm -> Write();
   hSm_mtdL_dLV_Lm -> Write();
   hLm_dLV -> Write();
   hSm_dLV -> Write();
   hLm_dLV_Lm -> Write();
   hSm_dLV_Lm -> Write();
   hSm_mtdS -> Write();
   hSm_mtdS_pvz -> Write();

   cLsbg -> Write();
   cSsbg -> Write();
   
   hMMABC -> Write();
   hMMABCpip -> Write();
   hMMABCpim -> Write();
   hinvMAC -> Write();
   hinvMAB -> Write();
   hMMAB -> Write();
   hMMBC -> Write();
   hDN -> Write();

   heventNo -> Write();
     
   fout -> Close();
}
