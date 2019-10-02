#include <TChain.h>
#include "TFile.h"
#include "TH1F.h"
#include <TLine.h>
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
#include "peakBG.C"

void nice_canv1(TVirtualPad * c)
{
    c -> SetBottomMargin(0.1);
    c -> SetLeftMargin(0.1);
    gStyle -> SetTitleSize(0.04, "xy");
    gStyle -> SetLabelSize(0.04, "xy");
//    gStyle -> SetTitleOffset(1.2, "x");
//    gStyle -> SetTitleOffset(1.3, "y");
}

void sigmaAna(){
    double LAMBDA_MASS = 1116; //MeV
    double SIGMA_MASS = 1382.8; //MeV
    
    //def histos
    int lm1 = 1060;
    int lm2 = 1300;
    int sm1 = 1200;
    int sm2 = 1800;
    int nbinsm = 200;
    int nbins3 = 100;
    int mtdmin = -5;
    int mtdmax = 100;
    
   TH1F *hLm = new TH1F("hLm", "#Lambda mass", nbinsm, lm1, lm2);
   TH1F *hpm = new TH1F("hpm", "proton mass", 200, 935, 940);
   TH1F *hpimm = new TH1F("hpimm", "#pion^{-} mass", 200, 130, 150);
   TH1F *hpipm = new TH1F("hpipm", "#pion^{+} mass", 200, 130, 150);
   TH1F *hSm = new TH1F("hSm", "#Sigma mass", nbinsm, sm1, sm2);

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

   TH1F *hMinTrackDist = new TH1F("hMinTrackDist", "MinTrackDist p-#pi^{-}", nbins3, mtdmin, mtdmax);
   TH1F *hMinTrackDistLpi = new TH1F("hMinTrackDistLpi", "MinTrackDist #Lambda-#pi^{+}", nbins3, mtdmin, mtdmax);
   TH2F *hMTD2D = new TH2F("hMTD2D", "MTD_{p-#pi^{-}} vs. MTD_{#Lambda-#pi^{+}}", nbins3, mtdmin, mtdmax, nbins3, mtdmin, mtdmax);
   TH1F *hPrimVertR = new TH1F("hPrimVertR", "Primary vertex R", 40, -10, 30);
   TH1F *hLDecVertR = new TH1F("hLDecVertR", "Lambda decay vertex R", 40, -10, 30);
   TH1F *hdVert = new TH1F("hdVert", "Distance between verticies (no cuts)", 100, -10, 90);
   TH1F *hdLVR = new TH1F("hdLVR", "Lambda decay vertex - Primary Vertex", 40, -10, 30);
   TH1F *hPrimVertZ = new TH1F("hPrimVertZ", "Primary vertex Z", 100, -100, 150);
   TH1F *hLDecVertX = new TH1F("hLDecVertX", "Lambda decay vertex X", 100, -100, 150);
   TH1F *hLDecVertY = new TH1F("hLDecVertY", "Lambda decay vertex Y", 100, -100, 150);
   TH1F *hLDecVertZ = new TH1F("hLDecVertZ", "Lambda decay vertex Z", 100, -100, 150);

   //mtd L cut
   TH1F *hLm_mtdL = new TH1F("hLm_mtdL", "#Lambda mass after MTD p-#pi^{-}", nbinsm, lm1, lm2);
   TH1F *hpipm_mtdL = new TH1F("hpipm_mtdL", "#pion^{+} mass after MTD p-#pi^{-}", 200, 130, 150);
   TH1F *hSm_mtdL = new TH1F("hSm_mtdL", "#Sigma mass after MTD p-#pi^{-}", nbinsm, sm1, sm2);

   //
   //mtd L & dVert cut
   TH1F *hLm_mtdL_dvert = new TH1F("hLm_mtdL_dvert", "Lambda mass after MTD p-#pi^{-} and dVert", nbinsm, lm1, lm2);
   TH1F *hSm_mtdL_dvert = new TH1F("hSm_mtdL_dvert", "Sigma mass after MTD #Lambda-#pi^{+} and dVert", nbinsm, sm1, sm2);
   //mtd L & dVert cut &Lm
   TH1F *hLm_mtdL_dvert_Lm = new TH1F("hLm_mtdL_dvert_Lm", "Lambda mass after MTD p-#pi^{-}, dVert and Lmass", nbinsm, lm1, lm2);
   TH1F *hSm_mtdL_dvert_Lm = new TH1F("hSm_mtdL_dvert_Lm", "Sigma mass after MTD #Lambda-#pi^{+}, dVert and Lmass", nbinsm, sm1, sm2);

   //mtd L & dVert cut & dLV cut
   TH1F *hLm_mtdL_dvert_dLV = new TH1F("hLm_mtdL_dvert_dLV", "Lambda mass after MTD p-#pi^{-}, dVert and (LDecVertR-PrimVertR)", nbinsm, lm1, lm2);
   TH1F *hpipm_mtdL_dvert_dLV = new TH1F("hpipm_mtdL_dvert_dLV", "pion+ mass after MTD p-#pi^{-}, dVert and (LDecVertR-PrimVertR)", nbinsm, lm1, lm2);
   TH1F *hSm_mtdL_dvert_dLV = new TH1F("hSm_mtdL_dvert_dLV", "Sigma mass after MTD #Lambda-#pi^{+}, dVert and (LDecVertR-PrimVertR)", nbinsm, sm1, sm2);

   //mtd L & dVert cut & dLV cut & Lm
   TH1F *hLm_mtdL_dvert_dLV_Lm = new TH1F("hLm_mtdL_dvert_dLV_Lm", "Lambda mass after MTD p-#pi^{-}, dVert, (LDecVertR-PrimVertR) and Lmass", nbinsm, lm1, lm2);
   TH1F *hSm_mtdL_dvert_dLV_Lm = new TH1F("hSm_mtdL_dvert_dLV_Lm", "Sigma mass after MTD #Lambda-#pi^{+}, dVert, (LDecVertR-PrimVertR) and Lmass", nbinsm, sm1, sm2);
   //
   //mtd L & dLV cut
   TH1F *hLm_mtdL_dLV = new TH1F("hLm_mtdL_dLV", "Lambda mass after MTD p-#pi^{-} and (LDecVertR-PrimVertR)", nbinsm, lm1, lm2);
   TH1F *hpipm_mtdL_dLV = new TH1F("hpipm_mtdL_dLV", "pion+ mass after MTD p-#pi^{-} and (LDecVertR-PrimVertR)", nbinsm, lm1, lm2);
   TH1F *hSm_mtdL_dLV = new TH1F("hSm_mtdL_dLV", "Sigma mass after MTD #Lambda-#pi^{+} and (LDecVertR-PrimVertR)", nbinsm, sm1, sm2);

   //mtd L & dLV cut & Lm
   TH1F *hLm_mtdL_dLV_Lm = new TH1F("hLm_mtdL_dLV_Lm", "Lambda mass after MTD p-#pi^{-}, (LDecVertR-PrimVertR) and Lmass", nbinsm, lm1, lm2);
   TH1F *hSm_mtdL_dLV_Lm = new TH1F("hSm_mtdL_dLV_Lm", "Sigma mass after MTD #Lambda-#pi^{+}, (LDecVertR-PrimVertR) and Lmass", nbinsm, sm1, sm2);
   //
   //dLV cut
   TH1F *hLm_dLV = new TH1F("hLm_dLV", "Lambda mass after (LDecVertR-PrimVertR)", nbinsm, lm1, lm2);
   TH1F *hpipm_dLV = new TH1F("hpipm_dLV", "pion+ mass after (LDecVertR-PrimVertR)", nbinsm, lm1, lm2);
   TH1F *hSm_dLV = new TH1F("hSm_dLV", "Sigma mass after (LDecVertR-PrimVertR)", nbinsm, sm1, sm2);

   //dLV cut & Lm
   TH1F *hLm_dLV_Lm = new TH1F("hLm_dLV_Lm", "Lambda mass after (LDecVertR-PrimVertR) and Lmass", nbinsm, lm1, lm2);
   TH1F *hSm_dLV_Lm = new TH1F("hSm_dLV_Lm", "Sigma mass after (LDecVertR-PrimVertR) and Lmass", nbinsm, sm1, sm2);

   //S: all L cuts &:
   //mtdS
   TH1F *hSm_mtdS = new TH1F("hSm_mtdS", "#Sigma mass after all cuts on #Lambda and MTD #Lambda-#pi^{+}", nbinsm, sm1, sm2);
   //mtdS & primVertZ
   TH1F *hSm_mtdS_pvz = new TH1F("hSm_mtdS_pvz", "#Sigma mass after all cuts on #Lambda, MTD #Lambda-#pi^{+} and PrimVertZ", nbinsm, sm1, sm2);
   
   //histos BG
   TH1F *hLmRecBG = new TH1F("hLmRecBG", "hLmRecBG", nbinsm, lm1, lm2);
   TH1F *hLmRecBGMtd = new TH1F("hLmRecBGMtd", "hLmRecBGMtd", nbinsm, lm1, lm2);
   TH1F *hLmRecBGMtdVertLz = new TH1F("hLmRecBGMtdVertLz", "hLmRecBGMtdVertLz", nbinsm, lm1, lm2);
   TH1F *hLmRecBGMtdDlv = new TH1F("hLmRecBGMtdDlv", "hLmRecBGMtdDlv", nbinsm, lm1, lm2);
   TH1F *hLmRecBGDlv = new TH1F("hLmRecBGDlv", "hLmRecBGDlv", nbinsm, lm1, lm2);
//   int cntLmRec, cntLmRecBG, cntLmRecMtd, cntLmRecMtdBG, cntLmRecMtdVertLz, cntLmRecMtdVertLzBG, cntLmRecMtdDlv, cntLmRecMtdDlvBG, cntLmRecDlv, cntLmRecDlvBG;
   //canvases with S&BG
   TCanvas *cLsbg1 = new TCanvas("cLsbg1", "cLsbg1", 1200, 800);
   TCanvas *cLsbg2 = new TCanvas("cLsbg2", "cLsbg2", 1200, 800);
   TCanvas *cSsbg1 = new TCanvas("cSsbg1", "cSsbg1", 1200, 800);
   TCanvas *cSsbg2 = new TCanvas("cSsbg2", "cSsbg2", 1200, 800);
   cLsbg1 -> Divide(2,1);
   cLsbg2 -> Divide(2,1);
   cSsbg1 -> Divide(2,1);
   cSsbg2 -> Divide(3,1);
   
   //miss mass spectra
   TH1F *hMMABC = new TH1F("hMMABC", "hMMABC", 100, -2300, 2350);
   TH1F *hMMABCpip = new TH1F("hMMABCpip", "hMMABCpip", 100, -3600, 2400);
   TH1F *hMMABCpim = new TH1F("hMMABCpim", "hMMABCpim", 100, -3600, 2400);
   TH1F *hMMAB = new TH1F("hMMAB", "hMMAB", 100, -2300, 2350);
   TH1F *hMMBC = new TH1F("hMMBC", "hMMBC", 100, -2300, 2350);
   //inv mass spectra
   TH1F *hinvMAC = new TH1F("hinvMAC", "hinvMAC", 200, 1055, 4000);
   TH1F *hinvMAB = new TH1F("hinvMAB", "hinvMAB", 200, 1066, 4000);

   int dminx = -2300;
   int dmaxx = 2350;
   int dminy = 1000;
   int dmaxy = 2500;
   TH2F *hDN = new TH2F("hDN", "hDN", 100, 1066, 1212, 100, -3600, 2400);
   TH2F *hDpim = new TH2F("hDpim", "hDpim", nbins3, dminx, dmaxx, nbins3, dminy, dmaxy);
   TH2F *hDpip = new TH2F("hDpip", "hDpip", nbins3, dminx, dmaxx, nbins3, dminy, dmaxy);
   TCanvas *cDpip = new TCanvas("cDpip", "cDpip");
   
   
   TH1F *hLm0 = new TH1F("hLm0", "#Lambda mass (no cuts)", nbinsm, lm1, lm2);
   TH1F *hSm0 = new TH1F("hSm0", "#Sigma mass (no cuts)", nbinsm, sm1, sm2);
   TH1F *hMinTrackDist0 = new TH1F("hMinTrackDist0", "MinTrackDist p-#pi^{-} (no cuts)", nbins3, mtdmin, mtdmax);
   TH1F *hMinTrackDistLpi0 = new TH1F("hMinTrackDistLpi0", "MinTrackDist #Lambda-#pi^{+} (no cuts)", nbins3, mtdmin, mtdmax);
   TH2F *hMTD2D0 = new TH2F("hMTD2D0", "MTD_{p-#pi^{-}} vs. MTD_{#Lambda-#pi^{+}} (no cuts)", nbins3, mtdmin, mtdmax, nbins3, mtdmin, mtdmax);
   TH1F *hPrimVertR0 = new TH1F("hPrimVertR0", "Primary vertex R (no cuts)", 40, -10, 30);
   TH1F *hLDecVertR0 = new TH1F("hLDecVertR0", "Lambda decay vertex R (no cuts)", 40, -10, 30);
   TH1F *hdLVR0 = new TH1F("hdLVR0", "Lambda decay vertex - Primary Vertex (no cuts)", 40, -10, 30);
   TH1F *hPrimVertZ0 = new TH1F("hPrimVertZ0", "Primary vertex Z (no cuts)", 100, -100, 150);
   TH1F *hLDecVertZ0 = new TH1F("hLDecVertZ0", "Lambda decay vertex Z (no cuts)", 100, -100, 150);
   TH1F *hdVert0 = new TH1F("hdVert0", "Distance between verticies (no cuts)", 100, -10, 90);
   TH1F *hMMABC_0 = new TH1F("hMMABC_0", "hMMABC (no cuts, no pid)", 100, -2300, 2350);
   TH1F *hMMABC0 = new TH1F("hMMABC0", "hMMABC (no cuts)", 100, -2300, 2350);
   TH1F *hinvMAC0 = new TH1F("hinvMAC0", "hinvMAC (no cuts)", 200, 0, 2200);
   TH1F *hinvMAC2 = new TH1F("hinvMAC2", "hinvMAC2 (no cuts)", 200, 0, 2200);
   TH1F *hinvMAB0 = new TH1F("hinvMAB0", "hinvMAB (no cuts)", 200, 1066, 4000);
   TH2F *hDpim0 = new TH2F("hDpim0", "hDpim (no cuts)", nbins3, dminx, dmaxx, nbins3, dminy, dmaxy);
   TH2F *hDpip0 = new TH2F("hDpip0", "hDpip (no cuts)", nbins3, dminx, dmaxx, nbins3, dminy, dmaxy);

   //side band
   TH1F *hLm_sb = new TH1F("hLm_sb", "hLm_sb", nbinsm, lm1, lm2);
   TH1F *hSm_sb = new TH1F("hSm_sb", "hSm_sb", nbinsm, sm1, sm2);
   int sbl1 = 1085;
   int sbl2 = 1100;
   int sbp1 = 1130;
   int sbp2 = 1150;
   
   //vertex studies
   int dvertmin = -100;
   int dvertmax = 100;
   int vertmin = -300;
   int vertmax = 300;
   int nbins1 = 200;
   int nbins2 = 400;
   int mommin = -1000;
   int mommax = 1000;
   int momzmin = 0;
   int momzmax = 3000;
   int dmommin = -400;
   int dmommax = 400;
   
/*   //vert_sim - vert_rec
   TH1F* hdvertlX = new TH1F("hdvertlX", "d(vertl_sim-vertl_rec)x", nbins2, dvertmin, dvertmax);
   TH1F* hdvertlY = new TH1F("hdvertlY", "d(vertl_sim-vertl_rec)y", nbins2, dvertmin, dvertmax);
   TH1F* hdvertlZ = new TH1F("hdvertlZ", "d(vertl_sim-vertl_rec)z", nbins2, dvertmin, dvertmax);
   TH1F* hdvertlR = new TH1F("hdvertlR", "d(vertl_sim-vertl_rec)r", nbins2, dvertmin, dvertmax);

   TH1F* hvx_lambdaX = new TH1F("hvx_lambdaX", "(vertl_rec)x", nbins2, vertmin, vertmax);
   TH1F* hvx_lambdaY = new TH1F("hvx_lambdaY", "(vertl_rec)y", nbins2, vertmin, vertmax);
   TH1F* hvx_lambdaZ = new TH1F("hvx_lambdaZ", "(vertl_rec)z", nbins2, vertmin, vertmax);

   TH1F* hgeantxvertexA = new TH1F("hgeantxvertexA", "(vertl_sim)x", nbins2, vertmin, vertmax);
   TH1F* hgeantyvertexA = new TH1F("hgeantyvertexA", "(vertl_sim)y", nbins2, vertmin, vertmax);
   TH1F* hgeantzvertexA = new TH1F("hgeantzvertexA", "(vertl_sim)z", nbins2, vertmin, vertmax);
   TH1F* hgeantxvertexB = new TH1F("hgeantxvertexB", "(vertl_sim)x", nbins2, vertmin, vertmax);
   TH1F* hgeantyvertexB = new TH1F("hgeantyvertexB", "(vertl_sim)y", nbins2, vertmin, vertmax);
   TH1F* hgeantzvertexB = new TH1F("hgeantzvertexB", "(vertl_sim)z", nbins2, vertmin, vertmax);
   
   TH1F* hvertex_recoGeant_x = new TH1F("hvertex_recoGeant_x", "vertex_recoGeant_x", nbins2, -1000, 1000);
   TH1F* hvertex_recoGeant_y = new TH1F("hvertex_recoGeant_y", "vertex_recoGeant_y", nbins2, dvertmin, dvertmax);
   TH1F* hvertex_recoGeant_z = new TH1F("hvertex_recoGeant_z", "vertex_recoGeant_z", nbins2, dvertmin, dvertmax);

   //part A
   TH1F* hMomAResx = new TH1F("hMomAResx", "MomAResx", nbins1, dmommin, dmommax);
   TH1F* hMomAResy = new TH1F("hMomAResy", "MomAResy", nbins1, dmommin, dmommax);
   TH1F* hMomAResz = new TH1F("hMomAResz", "MomAResz", nbins1, dmommin, dmommax);
   
   TH1F* hMomArecox = new TH1F("hMomArecox", "MomArecox", nbins1, mommin, mommax);
   TH1F* hMomArecoy = new TH1F("hMomArecoy", "MomArecoy", nbins1, mommin, mommax);
   TH1F* hMomArecoz = new TH1F("hMomArecoz", "MomArecoz", nbins1, momzmin, momzmax);

   TH1F* hMomAGeax = new TH1F("hMomAGeax", "MomAGeax", nbins1, mommin, mommax);
   TH1F* hMomAGeay = new TH1F("hMomAGeay", "MomAGeay", nbins1, mommin, mommax);
   TH1F* hMomAGeaz = new TH1F("hMomAGeaz", "MomAGeaz", nbins1, momzmin, momzmax);

   //part B
   TH1F* hMomBResx = new TH1F("hMomBResx", "MomBResx", nbins1, dmommin, dmommax);
   TH1F* hMomBResy = new TH1F("hMomBResy", "MomBResy", nbins1, dmommin, dmommax);
   TH1F* hMomBResz = new TH1F("hMomBResz", "MomBResz", nbins1, dmommin, dmommax);
   
   TH1F* hMomBrecox = new TH1F("hMomBrecox", "MomBrecox", nbins1, mommin, mommax);
   TH1F* hMomBrecoy = new TH1F("hMomBrecoy", "MomBrecoy", nbins1, mommin, mommax);
   TH1F* hMomBrecoz = new TH1F("hMomBrecoz", "MomBrecoz", nbins1, momzmin, momzmax);

   TH1F* hMomBGeax = new TH1F("hMomBGeax", "MomBGeax", nbins1, mommin, mommax);
   TH1F* hMomBGeay = new TH1F("hMomBGeay", "MomBGeay", nbins1, mommin, mommax);
   TH1F* hMomBGeaz = new TH1F("hMomBGeaz", "MomBGeaz", nbins1, momzmin, momzmax);

   //Lambda
   TH1F* hMomResx = new TH1F("hMomResx", "MomResx", nbins1, dmommin, dmommax);
   TH1F* hMomResy = new TH1F("hMomResy", "MomResy", nbins1, dmommin, dmommax);
   TH1F* hMomResz = new TH1F("hMomResz", "MomResz", nbins1, dmommin, dmommax);
   
   TH1F* hMomxLreco = new TH1F("hMomxLreco", "MomxLreco", nbins1, mommin, mommax);
   TH1F* hMomyLreco = new TH1F("hMomyLreco", "MomyLreco", nbins1, mommin, mommax);
   TH1F* hMomzLreco = new TH1F("hMomzLreco", "MomzLreco", nbins1, momzmin, momzmax);

   TH1F* hMomGeax = new TH1F("hMomGeax", "MomGeax", nbins1, mommin, mommax);
   TH1F* hMomGeay = new TH1F("hMomGeay", "MomGeay", nbins1, mommin, mommax);
   TH1F* hMomGeaz = new TH1F("hMomGeaz", "MomGeaz", nbins1, momzmin, momzmax);
*/

   TH1F* hrA = new TH1F("hrA", "hrA", nbins2, -100, 100);
   TH1F* hrB = new TH1F("hrB", "hrB", nbins2, -100, 100);

   TH1F* hLmass2 = new TH1F("hLmass2", "hLmass2", nbins2, 1000, 2000);
   TH1F* hrealLambdaM = new TH1F("hrealLambdaM", "hrealLambdaM", nbins2, 1085, 1140);
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
   float cmtd = 14;//25; //mm       
   float cvertLz1 = -20; //mm
   float cvertLz2 = 300; //mm
   float csVL1 = 15; //mm
   float csVL2 = 80; //mm
   float cdvert = 30; //mm distance between verticies
   float cLm1 = 1105; //MeV
   float cLm2 = 1125; //MeV
   float cmtdLpi = 15; //mm
   float cpvz1 = -55; //mm
   float cpvz2 = -30; //mm
//   float cmm1 = 1050; //miss mass ABC, MeV
   float cmm1 = 1350; //miss mass ABC, MeV
   float cmm2 = 1250; //inv mass AC, MeV
   
/*   //variables for the MTD L scan
   int mtdLi = 0;
   TH1F *hmtdLscan[15], *hmtdLscanSigma[15];
   TCanvas *cMtdScan = new TCanvas("cMtdScan", "cMtdScan", 1200, 800);
   cMtdScan -> Divide(5,3);
   TCanvas *cMtdScanSigma = new TCanvas("cMtdScanSigma", "cMtdScanSigma", 1200, 800);
   cMtdScanSigma -> Divide(5,3);
   char nameMtdLi[30], nameMtdLiSigma[30];
   double mtdVal[15], sgnVal[15];
   double sg, bg;
   for(int i = 0; i < 15; i++){
       mtdLi = (i+1)*2;
       sprintf(nameMtdLi, "hmtdLscan_%d", mtdLi);
       hmtdLscan[i] = new TH1F("hmtdLscan", nameMtdLi, 200, lm1, lm2);
       sprintf(nameMtdLiSigma, "hmtdLscanSigma_%d", mtdLi);
       hmtdLscanSigma[i] = new TH1F("hmtdLscanSigma", nameMtdLiSigma, 200, sm1, sm2);
   }
 
//variables for the MtdL(dVert) scan
   int dverti_mtdL = 0;
   int mtdLi_dvert = 0;
   TCanvas *cmtdLdvert[15], *cmtdLdvertSigma[15];
   TH1F *hdvertmtdLscan[15][15], *hdvertmtdLscanSigma[15][15];
   TGraph2D *gdVertmtdL = new TGraph2D();
   TH2F *hdVertmtdL = new TH2F("hdVertmtdL", "#alpha=S/#sqrt{S+B}", 15, 0, 75, 15, 0, 30);
   char namedverti_mtdL[30], namedverti_mtdLSigma[30];
   char nameCanv[30], nameCanvSigma[30];
   double sgnVal4[15][15];
   double dvertVal[15], sgnVal4a[15];
   double sg4, bg4;
   for(int i = 0; i < 15; i++){
       sprintf(nameCanv, "cmtdLdvert_%d", i);
       cmtdLdvert[i] = new TCanvas(nameCanv, nameCanv, 1200, 800);
       cmtdLdvert[i] -> Divide(5,3);
       sprintf(nameCanvSigma, "cmtdLdvertSigma_%d", i);
       cmtdLdvertSigma[i] = new TCanvas(nameCanvSigma, nameCanvSigma, 1200, 800);
       cmtdLdvertSigma[i] -> Divide(5,3);
       mtdLi_dvert = (i+1)*2;
       for(int j = 0; j < 15; j++){
	   dverti_mtdL = (j+1)*5;
	   sprintf(namedverti_mtdL, "hdvert_%d_mtdL_%d_scan", dverti_mtdL, mtdLi_dvert);
	   hdvertmtdLscan[i][j] = new TH1F("hdvertmtdLscan", namedverti_mtdL, 200, lm1, lm2);
	   sprintf(namedverti_mtdLSigma, "hdvert_%d_mtdL_%d_scanSigma", dverti_mtdL, mtdLi_dvert);
	   hdvertmtdLscanSigma[i][j] = new TH1F("hdvertmtdLscanSigma", namedverti_mtdL, 200, sm1, sm2);
       }
   }
*/ //<<<<end variables for the MTD L scan
   //variables for the MTD S scan                                                                                 
/*   int mtdLpii = 0;
   TH1F *hmtdLpiscan[15];
   TCanvas *cMtdLpiScan = new TCanvas("cMtdLpiScan", "cMtdLpiScan", 1200, 800);
   cMtdLpiScan -> Divide(5,3);
   char nameMtdLpii[30];
   double mtdLpiVal[15], sgnVal5[15];
   double sg5, bg5;
   for(int i = 0; i < 15; i++){
       mtdLpii = (i+1)*2;
       sprintf(nameMtdLpii, "hmtdLpiscan_%d", mtdLpii);
       hmtdLpiscan[i] = new TH1F("hmtdLpiscan", nameMtdLpii, 200, sm1, sm2);
   }                      
   */
   //read data
   TChain tree("T");
   //exp
   for(int i = 106; i < 128; i++){
      char inName[100];
       printf("apr07_day_%d.root\n",i);
       sprintf(inName, "/u/jkubos/analiza/gitdir/Sigma_JK/out_packed_3/apr07_day_%d.root", i);
//end exp
     
   //read no packed files
   //sim
   /*//  ifstream f("inList_all_001");
  // if(!f) cout << "FAILED TO OPEN DST FILE!" << endl;
  // const string line;
  // while(std::getline(f, line)){
   //    char inName[line.size() + 1];
  //     strcpy(inName, line.c_str());
   cout << "_" << endl ;
   int nSimCh = 2;
   for(int i = 1; i < nSimCh+1; i++){
       char inName[100];
       printf("output_00%d_all_0625.root\n",i);
       sprintf(inName, "/u/jkubos/analiza/gitdir/Sigma_JK/outputs_ch/output_00%d_all_0625.root", i);
*/     //end sim
       
     tree.Add(inName);
       
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
     primVertR, LDecVertR, primVertX, LDecVertX, primVertY, LDecVertY, primVertZ, LDecVertZ,
     mtd, mtdLpi,
       mmabc, mmabcpip, mmabcpim, invmac, invmab, mmab, mmbc,
       dvertlX, dvertlY, dvertlZ, dvertlR ,
       vx_lambdaX, vx_lambdaY, vx_lambdaZ,
       geantxvertexA, geantyvertexA, geantzvertexA,
       vertex_recogeant_x, vertex_recogeant_y, vertex_recogeant_z,
       momresx, momresy, momresz, momgeax, momgeay, momgeaz, momrecox, momrecoy, momrecoz,
       momaresx, momaresy, momaresz, momageax, momageay, momageaz, momarecox, momarecoy, momarecoz,
       mombresx, mombresy, mombresz, mombgeax, mombgeay, mombgeaz, mombrecox, mombrecoy, mombrecoz,
       r_hA, r_vB,
       realLambdaM, lmass2;
   Int_t geantIDa, geantIDb, geantIDc, geantParentIDa, geantParentIDb, geantParentIDc;
   

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
   tree.SetBranchAddress("fPrimaryVertexX", &primVertX);
   tree.SetBranchAddress("fPrimaryVertexY", &primVertY);
   tree.SetBranchAddress("fPrimaryVertexZ", &primVertZ);
   tree.SetBranchAddress("fLambdaDecayX", &LDecVertX);
   tree.SetBranchAddress("fLambdaDecayY", &LDecVertY);
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
   /* tree.SetBranchAddress("fdvertlX", &dvertlX);
   tree.SetBranchAddress("fdvertlY", &dvertlY);
   tree.SetBranchAddress("fdvertlZ", &dvertlZ);
   tree.SetBranchAddress("fDVres", &dvertlR);
   tree.SetBranchAddress("fvx_lambdaX", &vx_lambdaX);
   tree.SetBranchAddress("fvx_lambdaY", &vx_lambdaY);
   tree.SetBranchAddress("fvx_lambdaZ", &vx_lambdaZ);
   tree.SetBranchAddress("fgeantxvertexA", &geantxvertexA);
   tree.SetBranchAddress("fgeantyvertexA", &geantyvertexA);
   tree.SetBranchAddress("fgeantzvertexA", &geantzvertexA);
   tree.SetBranchAddress("fvertex_recoGeant_x", &vertex_recogeant_x);
   tree.SetBranchAddress("fvertex_recoGeant_y", &vertex_recogeant_y);
   tree.SetBranchAddress("fvertex_recoGeant_z", &vertex_recogeant_z);
   tree.SetBranchAddress("fMomAResx", &momaresx);
   tree.SetBranchAddress("fMomAResy", &momaresy);
   tree.SetBranchAddress("fMomAResz", &momaresz);
   tree.SetBranchAddress("fGeaAPx", &momageax);
   tree.SetBranchAddress("fGeaAPy", &momageay);
   tree.SetBranchAddress("fGeaAPz", &momageaz);
   tree.SetBranchAddress("fMomArecox", &momarecox);
   tree.SetBranchAddress("fMomArecoy", &momarecoy);
   tree.SetBranchAddress("fMomArecoz", &momarecoz);
   tree.SetBranchAddress("fMomResx", &momresx);
   tree.SetBranchAddress("fMomResy", &momresy);
   tree.SetBranchAddress("fMomResz", &momresz);
   tree.SetBranchAddress("fGeaPx", &momgeax);
   tree.SetBranchAddress("fGeaPy", &momgeay);
   tree.SetBranchAddress("fGeaPz", &momgeaz);
   tree.SetBranchAddress("fMomxLreco", &momrecox);
   tree.SetBranchAddress("fMomyLreco", &momrecoy);
   tree.SetBranchAddress("fMomzLreco", &momrecoz);
   tree.SetBranchAddress("fLmass2", &lmass2);*/
   tree.SetBranchAddress("fRA", &r_hA);
   tree.SetBranchAddress("fRB", &r_vB);
/*   tree.SetBranchAddress("frealLambdaM", &realLambdaM);
   tree.SetBranchAddress("fgpidA", &geantIDa);
   tree.SetBranchAddress("fgpidB", &geantIDb);
   tree.SetBranchAddress("fgpidC", &geantIDc);
   tree.SetBranchAddress("fgpidAparent", &geantParentIDa);
   tree.SetBranchAddress("fgpidBparent", &geantParentIDb);
   tree.SetBranchAddress("fgpidCparent", &geantParentIDc);
*/ 
   int eventNo = -1;
   TH1I *heventNo = new TH1I();
   heventNo -> SetName("heventNo");
   //event loop
   for (event = 0; event < nentries; event++){
     if(!(event%10000)){
       printf("Event: %ld\r",event);
     }
     tree.GetEntry(event); 

/*     if(Lm > 1105 && Lm < 1125){
	 //vert sim - vert reco
	 hdvertlX -> Fill(dvertlX);                                                                               
	 hdvertlY -> Fill(dvertlY);                                                                               
	 hdvertlZ -> Fill(dvertlZ);                                                                               
	 hdvertlR -> Fill(dvertlR);                                                                               
	 
	 hvx_lambdaX -> Fill(vx_lambdaX);
	 hvx_lambdaY -> Fill(vx_lambdaY);
	 hvx_lambdaZ -> Fill(vx_lambdaZ);
	 hgeantxvertexA -> Fill(geantxvertexA);
	 hgeantyvertexA -> Fill(geantyvertexA);
	 hgeantzvertexA -> Fill(geantzvertexA);
	 hvertex_recoGeant_x -> Fill(vertex_recogeant_x);
	 hvertex_recoGeant_y -> Fill(vertex_recogeant_y);
	 hvertex_recoGeant_z -> Fill(vertex_recogeant_z);
	 
	 hMomAResx -> Fill(momaresx);
	 hMomAResy -> Fill(momaresy);
	 hMomAResz -> Fill(momaresz);
	 hMomAGeax -> Fill(momageax);
	 hMomAGeay -> Fill(momageay);
	 hMomAGeaz -> Fill(momageaz);
	 hMomArecox -> Fill(momarecox);
	 hMomArecoy -> Fill(momarecoy);
	 hMomArecoz -> Fill(momarecoz);

	 hMomResx -> Fill(momresx);
	 hMomResy -> Fill(momresy);
	 hMomResz -> Fill(momresz);
	 hMomGeax -> Fill(momgeax);
	 hMomGeay -> Fill(momgeay);
	 hMomGeaz -> Fill(momgeaz);
	 hMomxLreco -> Fill(momrecox);
	 hMomyLreco -> Fill(momrecoy);
	 hMomzLreco -> Fill(momrecoz);
	 
	 hLmass2 -> Fill(lmass2);
	 
	 hrA -> Fill(r_hA);                                                                                       
	 hrB -> Fill(r_vB);   
	 
	 hrealLambdaM -> Fill(realLambdaM);
     }
*/
     float dVert = TMath::Sqrt( (LDecVertX - primVertX)*(LDecVertX - primVertX) + (LDecVertY - primVertY)*(LDecVertY - primVertY) + (LDecVertZ - primVertZ)*(LDecVertZ - primVertZ) );

//sim
     hMMABC_0 -> Fill(mmabc);
//     if(geantIDa != 14 || geantIDb != 9 || geantIDc != 8 || geantParentIDa != 18 || geantParentIDb != 18)
//	 continue;
     hLm0 -> Fill(Lm);
     hSm0 -> Fill(Sm);
     hPrimVertZ0 -> Fill(primVertZ);
     hLDecVertZ0 -> Fill(LDecVertZ);
     hPrimVertR0 -> Fill(primVertR);
     hLDecVertR0 -> Fill(LDecVertR);
     hdLVR0 -> Fill(LDecVertR - primVertR);
     hdVert0 -> Fill(dVert);
     hMinTrackDist0 -> Fill(mtd);
     hMinTrackDistLpi0 -> Fill(mtdLpi);
     hMTD2D0 -> Fill(mtd, mtdLpi);
     hMMABC0 -> Fill(mmabc);
     hinvMAC0 -> Fill(invmac);
     hinvMAB0 -> Fill(invmab);
     hDpim0 -> Fill(mmabc, invmab);
     hDpip0 -> Fill(mmabc, invmac);

     //cut on miss mass ABC
     if(mmabc < cmm1)
	 continue;
     hinvMAC2 -> Fill(invmac);
     
//     if(invmac > cmm2)
//	 continue;
   
     //no topological cuts
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
     hLDecVertX -> Fill(LDecVertX);
     hLDecVertY -> Fill(LDecVertY);
     hLDecVertZ -> Fill(LDecVertZ);
     hPrimVertR -> Fill(primVertR);
     hLDecVertR -> Fill(LDecVertR);
     hdLVR -> Fill(LDecVertR - primVertR);
     hdVert -> Fill(dVert);
     hMinTrackDist -> Fill(mtd);
     hMinTrackDistLpi -> Fill(mtdLpi);
     hMTD2D -> Fill(mtd, mtdLpi);
     hMMABC -> Fill(mmabc);
     hMMABCpip -> Fill(mmabcpip);
     hMMABCpim -> Fill(mmabcpim);
     hinvMAC -> Fill(invmac);
     hinvMAB -> Fill(invmab);
     hMMAB -> Fill(mmab);
     hMMBC -> Fill(mmbc);
     hDN -> Fill(invmac, mmabc);
     hDpim -> Fill(mmabc, invmab);
     hDpip -> Fill(mmabc, invmac);

/*     //MTD L && dvert scan
     for(int i = 0; i < 15; i++){
	 mtdLi = (i+1)*2;
	 if(mtd < mtdLi){
	     hmtdLscan[i] -> Fill(Lm);
	     if(Lm > cLm1 && Lm < cLm2)
		 hmtdLscanSigma[i] -> Fill(Sm);
	     //
	     for(int j = 0; j < 15; j++){
		 dverti_mtdL = (j+1)*5;
		 if(dVert > dverti_mtdL){
		     hdvertmtdLscan[i][j] -> Fill(Lm);
		     if(Lm > cLm1 && Lm < cLm2)
			 hdvertmtdLscanSigma[i][j] -> Fill(Sm);
		 }
	     }
	 }
     }
*/   
   //<<<<MTD L scan    
    //MTD S scan
     /*   if(mtd < cmtd && dVert > cdvert && Lm > cLm1 && Lm < cLm2){
	 for(int i = 0; i < 15; i++){
	     mtdLpii = (i+1)*2;
	     if(mtdLpi < mtdLpii)
		 hmtdLpiscan[i] -> Fill(Sm);
	 }
     }
     */
     //L & S inv mass histos
     //MTD L
     if(mtd < cmtd){
       hLm_mtdL -> Fill(Lm);
       hSm_mtdL -> Fill(Sm);
       //dVert
       if(dVert > cdvert){
	   hLm_mtdL_dvert -> Fill(Lm);
	   hSm_mtdL_dvert -> Fill(Sm);
        //LDecVertR-PrimVertR
/*       if((LDecVertR - primVertR) < csVL1){// && (LDecVertR - primVertR) < csVL2){
	   hLm_mtdL_dLV -> Fill(Lm);
	   hSm_mtdL_dLV -> Fill(Sm);
	   //vert sim - vert reco
	   hdvertlX -> Fill(dvertlX);
	   hdvertlY -> Fill(dvertlY);
	   hdvertlZ -> Fill(dvertlZ);
	   hdvertlR -> Fill(dvertlR);
	   
	   hrA -> Fill(r_hA);
	   hrB -> Fill(r_vB);
*/
	   if(Lm > cLm1 && Lm < cLm2){
	       hLm_mtdL_dvert_Lm -> Fill(Lm);
	       hSm_mtdL_dvert_Lm -> Fill(Sm);

	       if(mtdLpi < cmtdLpi){
		   hSm_mtdS -> Fill(Sm);
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

/*   //MTD L scan
   double partmpL[15];
   TF1 *ftmp, *fbgtmp, *fsigtmp;
   TH1F *htmp;
   char namePeaktmp[50];
   double a = 1085;
   double b = 1150;
   double fa = 1105;
   double fb = 1125;
   for(int i = 0; i < 15; i++){
       mtdLi = (i+1)*2;
       //calc Sgn for each mtdL value
       cMtdScan -> cd(i+1);
       char namehtmp[50];
       sprintf(namehtmp, "%htmp_mtdi%d", mtdLi);
       htmp = (TH1F*)hmtdLscan[i] -> Clone(namehtmp);
       const char *hnametmp = htmp -> GetName();
       namePeaktmp[] = 0;
       sprintf(namePeaktmp, "%s_peak", hnametmp);
       
       cout << "--" << i << "--" << endl;

       ftmp = new TF1("ftmp", "[0]*TMath::Voigt(x+[1],[2],[3]) + pol7(4)", a, b);
       fbgtmp = new TF1("fbgtmp", "pol7(0)", a, b);
       fsigtmp = new TF1("fsigtmp", "[0]*TMath::Voigt(x+[1],[2],[3])", a, b);

       if(i == 0){
	   ftmp -> SetParameters(
	       1.46482e+04,-1.11489e+03,1.35622e+00,1.61722e+00,
	       -4.21492e+06,3.65727e+03,3.32732e+00,6.71156e-05,-2.55953e-06,-2.31855e-09,2.01121e-12
	       );
       }else{
	   ftmp -> SetParameters(partmpL);
       }
       
       htmp -> Fit("ftmp", "0", "", a, b);
       htmp -> Fit("ftmp", "0", "", a, b);
       htmp -> SetName(namePeaktmp);

       partmpL[] = 0;
       ftmp -> GetParameters(partmpL);
       fsigtmp -> SetParameters(partmpL);
       fbgtmp -> SetParameters(&partmpL[4]);

       ftmp -> SetLineColor(kBlack);
       fsigtmp -> SetLineColor(kGreen);
       fbgtmp -> SetLineColor(kRed);
      
//       TH1F *histSigtmp = (TH1F*)htmp -> Clone("histSigtmp");
//       TH1F *histBGtmp = (TH1F*)htmp -> Clone("histBGtmp");
//       histSigtmp -> Add(fbgtmp, -1);
//       histSigtmp -> SetLineColor(kGreen);
//       histBGtmp -> Add(fsigtmp, -1);
//       histBGtmp -> SetLineColor(kRed);
//       histBGtmp -> SetMarkerColor(kRed);
//       histBGtmp -> SetMarkerStyle(20);
//       histBGtmp -> SetMarkerSize(.7);
       
//       double as = histSigtmp -> FindBin(1110);
//       double bs = histSigtmp -> FindBin(1120);
//       double ab = histBGtmp -> FindBin(1110);
//       double bb = histBGtmp -> FindBin(1120);
       sg = fsigtmp -> Integral(fa,fb);
       bg = fbgtmp -> Integral(fa,fb);
//     histSigtmp -> GetXaxis() -> SetRangeUser(1110,1120);
//     histBGtmp -> GetXaxis() -> SetRangeUser(1090,1140);
       
       htmp -> Draw();
       ftmp -> Draw("same");
       fbgtmp -> Draw("same");
       fsigtmp -> Draw("same");
       
       mtdVal[i] = mtdLi;
       sgnVal[i] = sg/sqrt(sg+bg);

       cMtdScanSigma -> cd(i+1);
       hmtdLscanSigma[i] -> Draw();
   }

   //MTD L_dVert scan
   int n = 0;
   double partmp[15];
   TF1 *ftmp, *fsigtmp, *fbgtmp;
   TH1F *htmp;
   char namePeaktmp[50];
   ofstream cout("output.txt");
   a = 1090;
   b = 1190;
   for(int i = 0; i < 15; i++){
       mtdLi_dvert = (i+1)*2;
       for(int j = 0; j < 15; j++){
	   cmtdLdvert[i] -> cd(j+1);
	   dverti_mtdL = (j+1)*5;
	   //calc Sgn for each mtdL value
	   char namehtmp[50];
	   sprintf(namehtmp, "%htmp_mtdi%d_dvert%d", mtdLi_dvert, dverti_mtdL);
	   htmp = (TH1F*)hdvertmtdLscan[i][j] -> Clone(namehtmp);
	   const char *hnametmp = htmp -> GetName();
	   namePeaktmp[] = 0;
	   sprintf(namePeaktmp, "%s_peak", hnametmp);
       
	   cout << "--mtdL" << i << "dvert" << j << "--" << endl;

	   ftmp = new TF1("ftmp", "[0]*TMath::Voigt(x+[1],[2],[3]) + pol6(4)", a, b);
	   fsigtmp = new TF1("fsigtmp", "[0]*TMath::Voigt(x+[1],[2],[3])", a, b);
	   fbgtmp = new TF1("fbgtmp", "pol6(0)", a, b);

	   double par0 = htmp -> GetBinContent(htmp -> GetMaximumBin());
	   if(j == 0){
	   ftmp -> SetParameters(
	       par0,-1.11489e+03,1.35622e+00,1.61722e+00,
	       -4.21492e+06,3.65727e+03,3.32732e+00,6.71156e-05,-2.55953e-06,-2.31855e-09,2.01121e-12
	       );
	   }else{
	       ftmp -> SetParameters(partmp);
	   }

	   
	   
	   htmp -> Fit("ftmp", "0", "", a, b);
	   htmp -> Fit("ftmp", "0", "", a, b);
	   htmp -> SetName(namePeaktmp);
	   
	   partmp[] = 0;
	   ftmp -> GetParameters(partmp);
	   fsigtmp -> SetParameters(partmp);
	   fbgtmp -> SetParameters(&partmp[4]);

//	   TH1F *histSigtmp = (TH1F*)htmp -> Clone("histSigtmp");
//	   TH1F *histBGtmp = (TH1F*)htmp -> Clone("histBGtmp");
//	   histSigtmp -> Add(fbgtmp, -1);
//	   histSigtmp -> SetLineColor(kGreen);
//	   histBGtmp -> Add(fsigtmp, -1);
//	   histBGtmp -> SetLineColor(kRed);
//	   histBGtmp -> SetMarkerColor(kRed);
//	   histBGtmp -> SetMarkerStyle(20);
//	   histBGtmp -> SetMarkerSize(.7);

	   ftmp -> SetLineColor(kBlack);
	   fsigtmp -> SetLineColor(kGreen);
	   fbgtmp -> SetLineColor(kRed);
	   
	   //double as = histSigtmp -> FindBin(1110);
	   //double bs = histSigtmp -> FindBin(1120);
	   //double ab = histBGtmp -> FindBin(1110);
	   //double bb = histBGtmp -> FindBin(1120);
	   sg4 = fsigtmp -> Integral(fa,fb);
	   bg4 = fbgtmp -> Integral(fa,fb);
	   //histSigtmp -> GetXaxis() -> SetRangeUser(1110,1120);
	   //histBGtmp -> GetXaxis() -> SetRangeUser(1090,1140);

	   htmp -> Draw();
	   ftmp -> Draw("same");
	   fsigtmp -> Draw("same");
	   fbgtmp -> Draw("same");

	   if(mtdLi_dvert == cmtd){
	       dvertVal[j] = dverti_mtdL;
	       if(sg4 < 0 || bg4 < 0)
		   sgnVal4a[j] = 0;
	       else
		   sgnVal4a[j] = sg4/sqrt(sg4+bg4);
	   }   
	       
	   if(sg4 < 0 || bg4 < 0)
	       sgnVal4[i][j] = 0;
	   else
	       sgnVal4[i][j] = sg4/sqrt(sg4+bg4);
	   cout << "sg4, bg4, n, mtdi, dvertj, sgn: " << sg4 << " " << bg4 << " " << n << " "<< dverti_mtdL << " " << mtdLi_dvert << " " << sgnVal4[i][j] << endl;
	   cout << ">>>>>>>>>>" << sgnVal4[i][j] << "<<<<<<<<<<<<<<" << endl;
	   for(int w = 0; w < 15; w++)
	       cout << partmp[w] << " " << endl;
	   
	   gdVertmtdL -> SetPoint(n, dverti_mtdL, mtdLi_dvert, sgnVal4[i][j]);
	   hdVertmtdL -> Fill(dverti_mtdL, mtdLi_dvert, sgnVal4[i][j]);
	   n++;

	   cmtdLdvertSigma[i] -> cd(j+1);
	   hdvertmtdLscanSigma[i][j] -> Draw();
       }
   }
*///<<<<<end scans
   //MTD S scan
       /*  double partmp[12];
   for(int i = 0; i < 15; i++){
       mtdLpii = (i+1)*2;
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
       partmp[] = 0;
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
       
     sg5 = histSigtmp -> Integral(fa,fb);
       bg5 = histBGtmp -> Integral(fa,fb);
       //    histSigtmp -> GetXaxis() -> SetRangeUser(1370,1405);
       //    histBGtmp -> GetXaxis() -> SetRangeUser(1370,1405);

       htmp -> Draw();
       histBGtmp -> Draw("same p");
       histSigtmp -> Draw("same");
       
       mtdLpiVal[i] = mtdLpii;
       sgnVal5[i] = sg5/sqrt(sg5+bg5);
       }*/

   //
/*   //<<<<   
   TGraph *gmtdscan = new TGraph(15, mtdVal, sgnVal);
   gmtdscan -> GetXaxis() -> SetTitle("MTD_{p#pi^{-}} [mm]");
   gmtdscan -> GetYaxis() -> SetTitle("#alpha = S/sqrt{S+B}");
   TGraph *gdVertscan = new TGraph(15, dvertVal, sgnVal4a);
   gdVertscan -> GetXaxis() -> SetTitle("dVert [mm]");
   gdVertscan -> GetYaxis() -> SetTitle("#alpha = S/sqrt{S+B}");

//   gdLVmtdL -> GetXaxis() -> SetTitle("LVert-PrimVert dist [mm]");
   gdVertmtdL -> GetXaxis() -> SetTitle("dVert [mm]");
   gdVertmtdL -> GetYaxis() -> SetTitle("MTD_{p#pi^{-}} [mm]");
   gdVertmtdL -> GetZaxis() -> SetTitle("#alpha = S/sqrt{S+B}");

//   TGraph *gmtdLpiscan = new TGraph(15, mtdLpiVal, sgnVal5);
//   gmtdLpiscan -> GetXaxis() -> SetTitle("MTD_{#Lambda#pi^{+}} [mm]");
//   gmtdLpiscan -> GetYaxis() -> SetTitle("#alpha = S/sqrt{S+B}");
 //<<<<<

   hdVertmtdL -> GetXaxis() -> SetTitle("dVert [mm]");
   hdVertmtdL -> GetYaxis() -> SetTitle("MTD_{p#pi^{-}} [mm]");
   hdVertmtdL -> GetZaxis() -> SetTitle("#alpha = S/sqrt{S+B}");
*/ 
   //BG subtruction
   //L no cuts
   double fita = 1100;
   double fitb = 1140;
   double fa = 1105;
   double fb = 1125;
   
   TH1F *histLnc0 = (TH1F*)hLm0 -> Clone("histLnc0");
   const char *hnameLnc0 = histLnc0 -> GetName();
   char namePeakLnc0[50];
   sprintf(namePeakLnc0, "%s_peak", hnameLnc0);

   cout << ">>>>>>>>fit BG no cuts<<<<<<<<<" << endl;
   TF1 * fhistLnc0 = new TF1("fhistLnc0", "[0]*TMath::Voigt(x+[1],[2],[3]) + pol6(4)", fita, fitb);
   TF1 * fsigLnc0 = new TF1("fsigLnc0", "[0]*TMath::Voigt(x+[1],[2],[3])", fita, fitb);
   TF1 * fbgLnc0 = new TF1("fbgLnc0", "pol6(0)", fita, fitb);

   double par0 = histLnc0 -> GetBinContent(histLnc0 -> GetMaximumBin());
   fhistLnc0 -> SetParameters(
       par0,-1.11489e+03,1.35622e+00,1.61722e+00,
       -4.21492e+06,3.65727e+03,3.32732e+00,6.71156e-05,-2.55953e-06,-2.31855e-09,2.01121e-12
       );

   histLnc0 -> Fit("fhistLnc0", "0", "", fita, fitb);
   histLnc0 -> Fit("fhistLnc0", "0", "", fita, fitb);
   double parLnc0[12];
   fhistLnc0 -> GetParameters(parLnc0);
   fsigLnc0 -> SetParameters(parLnc0);
   fbgLnc0 -> SetParameters(&parLnc0[4]);

   TH1F *histSigLnc0 = (TH1F*)histLnc0 -> Clone("histSigLnc0");
   TH1F *histBGLnc0 = (TH1F*)histLnc0 -> Clone("histBGLnc0");
   histSigLnc0 -> Add(fbgLnc0, -1);
   histSigLnc0 -> SetLineColor(kGreen);
   histBGLnc0 -> Add(fsigLnc0, -1);
   histBGLnc0 -> SetLineColor(kRed);
   histBGLnc0 -> SetMarkerColor(kRed);
   histBGLnc0 -> SetMarkerStyle(20);
   histBGLnc0 -> SetMarkerSize(.7);

   fsigLnc0 -> SetLineColor(kBlack);
   fbgLnc0 -> SetLineColor(kBlue);
          
   double cntLmSigLnc0 = fsigLnc0 -> Integral(fa,fb);
   double cntLmBGLnc0 = fbgLnc0 -> Integral(fa,fb);

   double a = histSigLnc0 -> FindBin(1110);
   double b = histSigLnc0 -> FindBin(1120);
   histSigLnc0 -> GetXaxis() -> SetRangeUser(1110,1120);
   histBGLnc0 -> GetXaxis() -> SetRangeUser(1090,1140);

   //L miss mass cut
   TH1F *histLnc = (TH1F*)hLm -> Clone("histLnc");
   const char *hnameLnc = histLnc -> GetName();
   char namePeakLnc[50];
   sprintf(namePeakLnc, "%s_peak", hnameLnc);

   cout << ">>>>>>>>fit BG no cuts<<<<<<<<<" << endl;
   TF1 * fhistLnc = new TF1("fhistLnc", "[0]*TMath::Voigt(x+[1],[2],[3]) + pol6(4)", fita, fitb);
   fhistLnc -> SetParameters(
       1.46482e+04,-1.11489e+03,1.35622e+00,1.61722e+00,
       -4.21492e+06,3.65727e+03,3.32732e+00,6.71156e-05,-2.55953e-06,-2.31855e-09,2.01121e-12
       );

   histLnc -> Fit("fhistLnc", "0", "", fita, fitb);
   histLnc -> Fit("fhistLnc", "0", "", fita, fitb);
   histLnc -> SetName(namePeakLnc);
   TF1 * fsigLnc = new TF1("fsigLnc", "[0]*TMath::Voigt(x+[1],[2],[3])", fita, fitb);
   TF1 * fbgLnc = new TF1("fbgLnc", "pol6(0)", fita, fitb);
   double parLnc[12];
   fhistLnc -> GetParameters(parLnc);
   fsigLnc -> SetParameters(parLnc);
   fbgLnc -> SetParameters(&parLnc[4]);

   TH1F *histSigLnc = (TH1F*)histLnc -> Clone("histSigLnc");
   TH1F *histBGLnc = (TH1F*)histLnc -> Clone("histBGLnc");
   histSigLnc -> Add(fbgLnc, -1);
   histSigLnc -> SetLineColor(kGreen);
   histBGLnc -> Add(fsigLnc, -1);
   histBGLnc -> SetLineColor(kRed);
   histBGLnc -> SetMarkerColor(kRed);
   histBGLnc -> SetMarkerStyle(20);
   histBGLnc -> SetMarkerSize(.7);

   fsigLnc->SetLineColor(kBlack);
   fbgLnc->SetLineColor(kBlue);
   
   double a = histSigLnc -> FindBin(1110);
   double b = histSigLnc -> FindBin(1120);
   histSigLnc -> GetXaxis() -> SetRangeUser(1110,1120);
   histBGLnc -> GetXaxis() -> SetRangeUser(1090,1140);

   double cntLmSigLnc = fsigLnc -> Integral(fa, fb);
   double cntLmBGLnc = fbgLnc -> Integral(fa, fb);
   /////

   //L mtd
   TH1F *histLmtd = (TH1F*)hLm_mtdL -> Clone("histLmtd");
   const char *hnameLmtd = histLmtd -> GetName();
   char namePeakLmtd[50]; 
   sprintf(namePeakLmtd, "%s_peak", hnameLmtd);

   cout << ">>>>>>>>fit BG MtdL cut<<<<<<<<<" << endl;
   TF1 * fhistLmtd = new TF1("fhistLmtd", "[0]*TMath::Voigt(x+[1],[2],[3]) + pol6(4)", fita, fitb);
   fhistLmtd -> SetParameters(
       1.46482e+04,-1.11489e+03,1.35622e+00,1.61722e+00,
       -4.21492e+06,3.65727e+03,3.32732e+00,6.71156e-05,-2.55953e-06,-2.31855e-09,2.01121e-12
       );

   histLmtd -> Fit("fhistLmtd", "0", "", fita, fitb);
   histLmtd -> Fit("fhistLmtd", "0", "", fita, fitb);
   histLmtd -> SetName(namePeakLmtd);
   TF1 * fsigLmtd = new TF1("fsigLmtd", "[0]*TMath::Voigt(x+[1],[2],[3])", fita, fitb);
   TF1 * fbgLmtd = new TF1("fbgLmtd", "pol6(0)", fita, fitb);
   double parLmtd[12];
   fhistLmtd -> GetParameters(parLmtd);
   fsigLmtd -> SetParameters(parLmtd);
   fbgLmtd -> SetParameters(&parLmtd[4]);

   TH1F *histSigLmtd = (TH1F*)histLmtd -> Clone("histSigLmtd");
   TH1F *histBGLmtd = (TH1F*)histLmtd -> Clone("histBGLmtd");
   histSigLmtd -> Add(fbgLmtd, -1);
   histSigLmtd -> SetLineColor(kGreen);
   histBGLmtd -> Add(fsigLmtd, -1);
   histBGLmtd -> SetLineColor(kRed);
   histBGLmtd -> SetMarkerColor(kRed);
   histBGLmtd -> SetMarkerStyle(20);
   histBGLmtd -> SetMarkerSize(.7);

   fsigLmtd->SetLineColor(kBlack);
   fbgLmtd->SetLineColor(kBlue);

   double a = histSigLmtd -> FindBin(1110);
   double b = histSigLmtd -> FindBin(1120);
   histSigLmtd -> GetXaxis() -> SetRangeUser(1110,1120);
   histBGLmtd -> GetXaxis() -> SetRangeUser(1090,1140);

   double cntLmSigLmtd = fsigLmtd -> Integral(fa, fb);
   double cntLmBGLmtd = fbgLmtd -> Integral(fa, fb);
   /////

   //L mtd dVert
   TH1F *histLmtdDvert = (TH1F*)hLm_mtdL_dvert -> Clone("histLmtdDvert");
   const char *hnameLmtdDvert = histLmtdDvert -> GetName();
   char namePeakLmtdDvert[50];
   sprintf(namePeakLmtdDvert, "%s_peak", hnameLmtdDvert);

   cout << ">>>>>>>>fit BG MtdL dVert cuts<<<<<<<<<" << endl;

   if(i < 5){
       fita = 1105;
       fitb = 1150;
   }else{
       fita = 1100;
       fitb = 1140;
   }
   
   
   double par0 = histLmtdDvert -> GetBinContent(histLmtdDvert -> GetMaximumBin());
   TF1 * fhistLmtdDvert = new TF1("fhistLmtdDvert", "[0]*TMath::Voigt(x+[1],[2],[3]) + pol6(4)", fita, fitb);
   fhistLmtdDvert -> SetParameters(
       //1.46482e+04,
       par0,-1.11489e+03,1.35622e+00,1.61722e+00,
       -4.21492e+06,3.65727e+03,3.32732e+00,6.71156e-05,-2.55953e-06,-2.31855e-09,2.01121e-12
       );

   histLmtdDvert -> Fit("fhistLmtdDvert", "0", "", fita, fitb);
   histLmtdDvert -> Fit("fhistLmtdDvert", "0", "", fita, fitb);
   histLmtdDvert -> SetName(namePeakLmtdDvert);
   TF1 * fsigLmtdDvert = new TF1("fsigLmtdDvert", "[0]*TMath::Voigt(x+[1],[2],[3]) ", fita, fitb);
   TF1 * fbgLmtdDvert = new TF1("fbgLmtdDvert", "pol6(0)", fita, fitb);
   double parLmtdDvert[12];
   fhistLmtdDvert -> GetParameters(parLmtdDvert);
   fsigLmtdDvert -> SetParameters(parLmtdDvert);
   fbgLmtdDvert -> SetParameters(&parLmtdDvert[4]);

   TH1F *histSigLmtdDvert = (TH1F*)histLmtdDvert -> Clone("histSigLmtdDvert");
   TH1F *histBGLmtdDvert = (TH1F*)histLmtdDvert -> Clone("histBGLmtdDvert");
   histSigLmtdDvert -> Add(fbgLmtdDvert, -1);
   histSigLmtdDvert -> SetLineColor(kGreen);
   histBGLmtdDvert -> Add(fsigLmtdDvert, -1);
   histBGLmtdDvert -> SetLineColor(kRed);
   histBGLmtdDvert -> SetMarkerColor(kRed);
   histBGLmtdDvert -> SetMarkerStyle(20);
   histBGLmtdDvert -> SetMarkerSize(.7);

   fsigLmtdDvert->SetLineColor(kBlack);
   fbgLmtdDvert->SetLineColor(kBlue);

   double a = histSigLmtdDvert -> FindBin(1110);
   double b = histSigLmtdDvert -> FindBin(1120);
   histSigLmtdDvert -> GetXaxis() -> SetRangeUser(1110,1120);
   histBGLmtdDvert -> GetXaxis() -> SetRangeUser(1090,1140);

   double cntLmSigLmtdDvert = fsigLmtdDvert -> Integral(fa, fb);
   double cntLmBGLmtdDvert = fbgLmtdDvert -> Integral(fa, fb);
   ///

   //>>>>>>>>>>>>>Sigma<<<<<<<<<<<<<<<<<<<<<<<
   //S no cuts
   int fitS1 = 1200;
   int fitS2 = 1700;
   int fitSsig1 = 1325;
   int fitSsig2 = 1445;
   
   TH1F *histSnc = (TH1F*)hSm -> Clone("histSnc");
   histSnc -> GetXaxis() -> SetTitle("m_{p#pi^{+}#pi^{-}}");
   histSnc -> GetXaxis() -> SetTitleSize(.05);
   histSnc -> GetXaxis() -> SetLabelSize(.05);
   histSnc -> Draw();
   
   const char *hnameSnc = histSnc -> GetName();
   char namePeakSnc[50];
   sprintf(namePeakSnc, "%s_peak", hnameSnc);

   cout << ">>>>>>>>fit BG no cuts<<<<<<<<<" << endl;
   TF1 * fhistSnc = new TF1("fhistSnc", "[0]*[1]*[1]/( (x*x-[2]*[2])*(x*x-[2]*[2]) + (x*x*x*x*[1]*[1]/([2]*[2])) ) + pol6(3)", fitS1, fitS2);
   fhistSnc -> SetParameters(
       1.40794e+09, 36, 1.38415e+03,
       -9.08964e+05,2.35372e+03,-2.33096e+00,
       1.08507e-03,-2.28651e-07,1.57016e-11
       );
   fhistSnc->SetParLimits(0,0,10000000000);
   fhistSnc->SetParLimits(1,30,40);
   fhistSnc->SetParLimits(2,1360,1400);
   
   histSnc -> Fit("fhistSnc", "0", "", fitS1, fitS2);
   histSnc -> Fit("fhistSnc", "0", "", fitS1, fitS2);
   histSnc -> SetName(namePeakSnc);
   TF1 * fsigSnc = new TF1("fsigSnc", "[0]*[1]*[1]/( (x*x-[2]*[2])*(x*x-[2]*[2]) + (x*x*x*x*[1]*[1]/([2]*[2])) )", fitSsig1, fitSsig2);
   TF1 * fbgSnc = new TF1("fbgSnc", "pol6(0)", fitS1, fitS2);
   double parSnc[12];
   fhistSnc -> GetParameters(parSnc);
   fsigSnc -> SetParameters(parSnc);
   fbgSnc -> SetParameters(&parSnc[3]);

   TH1F *histSigSnc = (TH1F*)histSnc -> Clone("histSigSnc");
   TH1F *histBGSnc = (TH1F*)histSnc -> Clone("histBGSnc");
   histSigSnc -> Add(fbgSnc, -1);
   histSigSnc -> SetLineColor(kGreen);
   histBGSnc -> Add(fsigSnc, -1);
   histBGSnc -> SetLineColor(kRed);
   histBGSnc -> SetMarkerColor(kRed);
   histBGSnc -> SetMarkerStyle(20);
   histBGSnc -> SetMarkerSize(.7);

   fsigSnc -> SetLineColor(kBlack);
   fsigSnc -> SetMarkerColor(kBlack);
   fsigSnc -> SetMarkerStyle(20);
   fsigSnc -> SetMarkerSize(.5);
   fsigSnc -> Draw("same");
   fbgSnc -> SetLineColor(kBlue);
   fbgSnc -> SetLineColor(kBlue);
   fbgSnc -> SetMarkerColor(kBlue);
   fbgSnc -> SetMarkerStyle(20);
   fbgSnc -> SetMarkerSize(.5);
   fbgSnc -> Draw("same");
   
   double cntSmSigSnc =  fsigSnc -> Integral(fitSsig1,fitSsig2);
   double cntSmBGSnc = fbgSnc -> Integral(fitSsig1,fitSsig2);
   /////
 
   //S mtdL
   TH1F *histSmtdL = (TH1F*)hSm_mtdL -> Clone("histSmtdL");
   histSmtdL -> GetXaxis() -> SetTitle("m_{p#pi^{+}#pi^{-}}");
   histSmtdL -> GetXaxis() -> SetTitleSize(.05);
   histSmtdL -> GetXaxis() -> SetLabelSize(.05);
   histSmtdL -> Draw();

   const char *hnameSmtdL = histSmtdL -> GetName();
   char namePeakSmtdL[50];
   sprintf(namePeakSmtdL, "%s_peak", hnameSmtdL);

   cout << ">>>>>>>>fit BG MtdL cut<<<<<<<<<" << endl;
   TF1 * fhistSmtdL = new TF1("fhistSmtdL", "[0]*[1]*[1]/( (x*x-[2]*[2])*(x*x-[2]*[2]) + (x*x*x*x*[1]*[1]/([2]*[2])) ) + pol6(3)", fitS1, fitS2);
   fhistSmtdL -> SetParameters(
       1.40794e+09, 36, 1.38415e+03,
       -9.08964e+05,2.35372e+03,-2.33096e+00,
       1.08507e-03,-2.28651e-07,1.57016e-11
       );

   fhistSmtdL -> SetParLimits(0,0,10000000000);
   fhistSmtdL -> SetParLimits(1,30,40);
   fhistSmtdL -> SetParLimits(2,1360,1400);
   
   histSmtdL -> Fit("fhistSmtdL", "0", "", fitS1, fitS2);
   histSmtdL -> Fit("fhistSmtdL", "0", "", fitS1, fitS2);
   histSmtdL -> SetName(namePeakSmtdL);
   TF1 * fsigSmtdL = new TF1("fsigSmtdL", "[0]*[1]*[1]/( (x*x-[2]*[2])*(x*x-[2]*[2]) + (x*x*x*x*[1]*[1]/([2]*[2])) )", fitSsig1,fitSsig2);
   TF1 * fbgSmtdL = new TF1("fbgSmtdL", "pol6(0)", fitS1, fitS2);
   double parSmtdL[12];
   fhistSmtdL -> GetParameters(parSmtdL);
   fsigSmtdL -> SetParameters(parSmtdL);
   fbgSmtdL -> SetParameters(&parSmtdL[3]);

   TH1F *histSigSmtdL = (TH1F*)histSmtdL -> Clone("histSigSmtdL");
   TH1F *histBGSmtdL = (TH1F*)histSmtdL -> Clone("histBGSmtdL");
   histSigSmtdL -> Add(fbgSmtdL, -1);
   histSigSmtdL -> SetLineColor(kGreen);
   histBGSmtdL -> Add(fsigSmtdL, -1);
   histBGSmtdL -> SetLineColor(kRed);
   histBGSmtdL -> SetMarkerColor(kRed);
   histBGSmtdL -> SetMarkerStyle(20);
   histBGSmtdL -> SetMarkerSize(.7);

   fsigSmtdL -> SetLineColor(kBlack);
   fsigSmtdL -> SetMarkerColor(kBlack);
   fsigSmtdL -> SetMarkerStyle(20);
   fsigSmtdL -> SetMarkerSize(.5);
   fsigSmtdL -> Draw("same");
   fbgSmtdL -> SetLineColor(kBlue);
   fbgSmtdL -> SetMarkerColor(kBlue);
   fbgSmtdL -> SetMarkerStyle(20);
   fbgSmtdL -> SetMarkerSize(.5);
   fbgSmtdL -> Draw("same");

   double cntSmSigSmtdL = fsigSmtdL -> Integral(fitSsig1,fitSsig2);
   double cntSmBGSmtdL = fbgSmtdL -> Integral(fitSsig1,fitSsig2);
   /////
 
  //S mtdL dVert
   TH1F *histSmtdLDvert = (TH1F*)hSm_mtdL_dvert -> Clone("histSmtdLDvert");
   histSmtdLDvert -> GetXaxis() -> SetTitle("m_{p#pi^{+}#pi^{-}}");
   histSmtdLDvert -> GetXaxis() -> SetTitleSize(.05);
   histSmtdLDvert -> GetXaxis() -> SetLabelSize(.05);
   histSmtdLDvert -> Draw();

   const char *hnameSmtdLDvert = histSmtdLDvert -> GetName();
   char namePeakSmtdLDvert[50];
   sprintf(namePeakSmtdLDvert, "%s_peak", hnameSmtdLDvert);
 
   cout << ">>>>>>>>fit BG MtdL dVert cuts<<<<<<<<<" << endl;
   TF1 * fhistSmtdLDvert = new TF1("fhistSmtdLDvert", "[0]*[1]*[1]/( (x*x-[2]*[2])*(x*x-[2]*[2]) + (x*x*x*x*[1]*[1]/([2]*[2])) ) + pol6(3)", fitS1, fitS2);
   fhistSmtdLDvert -> SetParameters(
       1.40794e+09, 36, 1.38415e+03,
       -9.08964e+05,2.35372e+03,-2.33096e+00,
       1.08507e-03,-2.28651e-07,1.57016e-11
       );
   fhistSmtdLDvert -> SetParLimits(0,0,10000000000);
   fhistSmtdLDvert -> SetParLimits(1,30,40);
   fhistSmtdLDvert -> SetParLimits(2,1360,1400);
   
   histSmtdLDvert -> Fit("fhistSmtdLDvert", "0", "", fitS1, fitS2);
   histSmtdLDvert -> Fit("fhistSmtdLDvert", "0", "", fitS1, fitS2);
   histSmtdLDvert -> SetName(namePeakSmtdLDvert);
   TF1 * fsigSmtdLDvert = new TF1("fsigSmtdLDvert", "[0]*[1]*[1]/( (x*x-[2]*[2])*(x*x-[2]*[2]) + (x*x*x*x*[1]*[1]/([2]*[2])) )", fitS1,fitS2);
   TF1 * fbgSmtdLDvert = new TF1("fbgSmtdLDvert", "pol6(0)", fitS1, fitS2);
   double parSmtdLDvert[12];
   fhistSmtdLDvert -> GetParameters(parSmtdLDvert);
   fsigSmtdLDvert -> SetParameters(parSmtdLDvert);
   fbgSmtdLDvert -> SetParameters(&parSmtdLDvert[3]);

   TH1F *histSigSmtdLDvert = (TH1F*)histSmtdLDvert -> Clone("histSigSmtdLDvert");
   TH1F *histBGSmtdLDvert = (TH1F*)histSmtdLDvert -> Clone("histBGSmtdLDvert");
   histSigSmtdLDvert -> Add(fbgSmtdLDvert, -1);
   histSigSmtdLDvert -> SetLineColor(kGreen);
   histBGSmtdLDvert -> Add(fsigSmtdLDvert, -1);
   histBGSmtdLDvert -> SetLineColor(kRed);
   histBGSmtdLDvert -> SetMarkerColor(kRed);
   histBGSmtdLDvert -> SetMarkerStyle(20);
   histBGSmtdLDvert -> SetMarkerSize(.7);

   fsigSmtdLDvert -> SetLineColor(kBlack);
   fsigSmtdLDvert -> SetMarkerColor(kBlack);
   fsigSmtdLDvert -> SetMarkerStyle(20);
   fsigSmtdLDvert -> SetMarkerSize(.5);
   fsigSmtdLDvert -> Draw("same");
   fbgSmtdLDvert -> SetLineColor(kBlue);
   fbgSmtdLDvert -> SetMarkerColor(kBlue);
   fbgSmtdLDvert -> SetMarkerStyle(20);
   fbgSmtdLDvert -> SetMarkerSize(.5);
   fbgSmtdLDvert -> Draw("same");
          
   double cntSmSigSmtdLDvert = fsigSmtdLDvert -> Integral(fitSsig1,fitSsig2);
   double cntSmBGSmtdLDvert = fbgSmtdLDvert -> Integral(fitSsig1,fitSsig2);
   ///

   //S mtdL dVert Lm
   //BW+pol7 fit
   TCanvas *cSigmaFit = new TCanvas("cSigmaFit", "cSigmaFit");
   cSigmaFit -> Divide(3,2);
   cSigmaFit -> cd(1);
   nice_canv1(gPad);
   hLm0 -> GetXaxis() -> SetTitleSize(.05);
   hLm0 -> GetXaxis() -> SetLabelSize(.05);
   hLm0 -> Draw();
   cSigmaFit -> cd(4);
   nice_canv1(gPad);
   hLm0 -> GetXaxis() -> SetTitleSize(.05);
   hLm0 -> GetXaxis() -> SetLabelSize(.05);
   hSm0 -> Draw();
   cSigmaFit -> cd(2);
   nice_canv1(gPad);
   hLm -> GetXaxis() -> SetTitleSize(.05);
   hLm -> GetXaxis() -> SetLabelSize(.05);
   hLm -> Draw();
   cSigmaFit -> cd(5);
   nice_canv1(gPad);
   hSm -> GetXaxis() -> SetTitleSize(.05);
   hSm -> GetXaxis() -> SetLabelSize(.05);
   hSm -> Draw();
   cSigmaFit -> cd(3);
   nice_canv1(gPad);
   hLm_mtdL_dLV -> GetXaxis() -> SetTitleSize(.05);
   hLm_mtdL_dLV -> GetXaxis() -> SetLabelSize(.05);
   hLm_mtdL_dLV -> Draw();
   cSigmaFit -> cd(6);

   TH1F *histSmtdLDvertLm = (TH1F*)hSm_mtdL_dvert_Lm -> Clone("histSmtdLDvertLm");
   histSmtdLDvertLm -> GetXaxis() -> SetTitle("m_{p#pi^{+}#pi^{-}}");
   histSmtdLDvertLm -> GetXaxis() -> SetTitleSize(.05);
   histSmtdLDvertLm -> GetXaxis() -> SetLabelSize(.05);
   histSmtdLDvertLm -> Draw();
   
   const char *hnameSmtdLDvertLm = histSmtdLDvertLm -> GetName();
   char namePeakSmtdLDvertLm[50];
   sprintf(namePeakSmtdLDvertLm, "%s_peak", hnameSmtdLDvertLm);

   cout << ">>>>>>>>fit BG MtdL dVert Lm cuts<<<<<<<<<" << endl;
   histSmtdLDvertLm -> SetName(namePeakSmtdLDvertLm);
   TF1 * fsigSmtdLDvertLm = new TF1("fsigSmtdLDvertLm", "[0]*[1]*[1]/( (x*x-[2]*[2])*(x*x-[2]*[2]) + (x*x*x*x*[1]*[1]/([2]*[2])) )", fitSsig1, fitSsig2);
   TF1 * fbgSmtdLDvertLm = new TF1("fbgSmtdLDvertLm", "pol6(0)", fitS1, fitS2);
   double parSmtdLDvertLm[12] = peakBG(histSmtdLDvertLm);
   fhistSmtdLDvertLm -> GetParameters(parSmtdLDvertLm);
   fsigSmtdLDvertLm -> SetParameters(parSmtdLDvertLm);
   fbgSmtdLDvertLm -> SetParameters(&parSmtdLDvertLm[3]);

   TH1F *histSigSmtdLDvertLm = (TH1F*)histSmtdLDvertLm -> Clone("histSigSmtdLDvertLm");
   TH1F *histBGSmtdLDvertLm = (TH1F*)histSmtdLDvertLm -> Clone("histBGSmtdLDvertLm");
   histSigSmtdLDvertLm -> Add(fbgSmtdLDvertLm, -1);
   histSigSmtdLDvertLm -> SetLineColor(kGreen);
   histBGSmtdLDvertLm -> Add(fsigSmtdLDvertLm, -1);
   histBGSmtdLDvertLm -> SetLineColor(kRed);
   histBGSmtdLDvertLm -> SetMarkerColor(kRed);
   histBGSmtdLDvertLm -> SetMarkerStyle(20);
   histBGSmtdLDvertLm -> SetMarkerSize(.7);

   fsigSmtdLDvertLm -> SetLineColor(kBlack);
   fsigSmtdLDvertLm -> SetMarkerColor(kBlack);
   fsigSmtdLDvertLm -> SetMarkerStyle(20);
   fsigSmtdLDvertLm -> SetMarkerSize(.5);
   fsigSmtdLDvertLm -> Draw("same");
   fbgSmtdLDvertLm -> SetLineColor(kBlue);
   fbgSmtdLDvertLm -> SetMarkerColor(kBlue);
   fbgSmtdLDvertLm -> SetMarkerStyle(20);
   fbgSmtdLDvertLm -> SetMarkerSize(.5);
   fbgSmtdLDvertLm -> Draw("same");

   TH1F *hpeakBG = (TH1F*)histSmtdLDvertLm -> Clone("hpeakBG");
   hpeakBG -> Add(fhistSmtdLDvertLm, -1);
   hpeakBG -> Draw("same");
      
   double cntSmSigSmtdLDvertLm = fsigSmtdLDvertLm -> Integral(fitSsig1,fitSsig2);
   double cntSmBGSmtdLDvertLm = fbgSmtdLDvertLm -> Integral(fitSsig1,fitSsig2);

   //end gaus+pol7 fit

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
   
   ///
   //S mtdL dVert Lm mtdS
   TH1F *histSmtdS = (TH1F*)hSm_mtdS -> Clone("histSmtdS");
   histSmtdS -> GetXaxis() -> SetTitle("m_{p#pi^{+}#pi^{-}}");
   histSmtdS -> GetXaxis() -> SetTitleSize(.05);
   histSmtdS -> GetXaxis() -> SetLabelSize(.05);
   histSmtdS -> Draw();
   
   const char *hnameSmtdS = histSmtdS -> GetName();
   char namePeakSmtdS[50];
   sprintf(namePeakSmtdS, "%s_peak", hnameSmtdS);
   
   cout << ">>>>>>>>fit BG MtdL dVert MtdS cuts<<<<<<<<<" << endl;
   TF1 * fhistSmtdS = new TF1("fhistSmtdS", "[0]*[1]*[1]/( (x*x-[2]*[2])*(x*x-[2]*[2]) + (x*x*x*x*[1]*[1]/([2]*[2])) ) + pol6(3)", fitS1, fitS2);
   fhistSmtdS -> SetParameters(
       1.40794e+09, 36, 1.38415e+03,
       -9.08964e+05,2.35372e+03,-2.33096e+00,
       1.08507e-03,-2.28651e-07,1.57016e-11
       );
   fhistSmtdS -> SetParLimits(0,0,10000000000);
   fhistSmtdS -> SetParLimits(1,30,40);
   fhistSmtdS -> SetParLimits(2,1360,1400);
   
   histSmtdS -> Fit("fhistSmtdS", "0", "", fitS1, fitS2);
   histSmtdS -> Fit("fhistSmtdS", "0", "", fitS1, fitS2);
   histSmtdS -> SetName(namePeakSmtdS);
   TF1 * fsigSmtdS = new TF1("fsigSmtdS", "[0]*[1]*[1]/( (x*x-[2]*[2])*(x*x-[2]*[2]) + (x*x*x*x*[1]*[1]/([2]*[2])) )", fitSsig1, fitSsig2);
   TF1 * fbgSmtdS = new TF1("fbgSmtdS", "pol6(0)", fitS1, fitS2);
   double parSmtdS[12];
   fhistSmtdS -> GetParameters(parSmtdS);
   fsigSmtdS -> SetParameters(parSmtdS);
   fbgSmtdS -> SetParameters(&parSmtdS[3]);

   TH1F *histSigSmtdS = (TH1F*)histSmtdS -> Clone("histSigSmtdS");
   TH1F *histBGSmtdS = (TH1F*)histSmtdS -> Clone("histBGSmtdS");
   histSigSmtdS -> Add(fbgSmtdS, -1);
   histSigSmtdS -> SetLineColor(kGreen);
   histBGSmtdS -> Add(fsigSmtdS, -1);
   histBGSmtdS -> SetLineColor(kRed);
   histBGSmtdS -> SetMarkerColor(kRed);
   histBGSmtdS -> SetMarkerStyle(20);
   histBGSmtdS -> SetMarkerSize(.7);

   fsigSmtdS -> SetLineColor(kBlack);
   fsigSmtdS -> SetMarkerColor(kBlack);
   fsigSmtdS -> SetMarkerStyle(20);
   fsigSmtdS -> SetMarkerSize(.5);
   fsigSmtdS -> Draw("same");
   fbgSmtdS -> SetLineColor(kBlue);
   fbgSmtdS -> SetMarkerColor(kBlue);
   fbgSmtdS -> SetMarkerStyle(20);
   fbgSmtdS -> SetMarkerSize(.5);
   fbgSmtdS -> Draw("same");
   
   double cntSmSigSmtdS = fsigSmtdS -> Integral(fitSsig1,fitSsig2);
   double cntSmBGSmtdS = fbgSmtdS -> Integral(fitSsig1,fitSsig2);

 ///
   //end BG
   
   cout << "nentries=" << nentries << endl << "S: " << endl << "Lambda:" << " Lnc=" << cntLmSigLnc << " LmtdL=" << cntLmSigLmtd << " LmtdLdVert=" << cntLmSigLmtdDvert << endl;
   cout << "Sigma: " << " Snc=" << cntSmSigSnc << " SmtdL=" << cntSmSigSmtdL << " SmtdLdVert=" << cntSmSigSmtdLDvert << " SLcuts=" << cntSmSigSmtdLDvertLm  << "SLcutsMtdS" << cntSmSigSmtdS << endl;
   
   cout << "B: " << endl << "Lambda:" << " Lnc=" << cntLmBGLnc << " LmtdL=" << cntLmBGLmtd << " LmtdLdVert=" << cntLmBGLmtdDvert << endl;
   cout << "Sigma: " << " Snc=" << cntSmBGSnc << " SmtdL=" << cntSmBGSmtdL << " SmtdLdVert=" << cntSmBGSmtdLDvert << " SLcuts=" << cntSmBGSmtdLDvertLm << " SLcutsMtdS" << cntSmBGSmtdS << endl;
    
    cout << "S/B: " << endl << "Lambda: " << " Lnc=" << cntLmSigLnc/cntLmBGLnc << " LmtdL=" << cntLmSigLmtd/cntLmBGLmtd << " LmtdLdVert=" << cntLmSigLmtdDvert/cntLmBGLmtdDvert << endl;
    cout << "Sigma: " << " Snc=" << cntSmSigSnc/cntSmBGSnc << " SmtdL=" << cntSmSigSmtdL/cntSmSigSmtdL << " SmtdLdVert=" << cntSmSigSmtdLDvert/cntSmBGSmtdLDvert << " SLcuts=" << cntSmSigSmtdLDvertLm/cntSmBGSmtdLDvertLm << " SLcutsMtdS" << cntSmSigSmtdS/cntSmBGSmtdS << endl;       

    double sgfLnc, sgfLmtdL, sgfLmtdLdVert,
	sgfSnc, sgfSmtdL, sgfSmtdLdVert, sgfSLcuts, sgfSLcutsMtdS;
    if((cntLmSigLnc + cntLmBGLnc) > 0)
	sgfLnc = cntLmSigLnc/sqrt(cntLmSigLnc + cntLmBGLnc);
    if((cntLmSigLmtd + cntLmBGLmtd) > 0)
	sgfLmtdL = cntLmSigLmtd/sqrt(cntLmSigLmtd + cntLmBGLmtd);
    if((cntLmSigLmtdDvert + cntLmBGLmtdDvert) > 0)
	sgfLmtdLdVert = cntLmSigLmtdDvert/sqrt(cntLmSigLmtdDvert + cntLmBGLmtdDvert);
    if((cntSmSigSnc + cntSmBGSnc) > 0)
	sgfSnc = cntSmSigSnc/sqrt(cntSmSigSnc + cntSmBGSnc);
    if((cntSmSigSmtdL + cntSmBGSmtdL) > 0)
	sgfSmtdL = cntSmSigSmtdL/sqrt(cntSmSigSmtdL + cntSmBGSmtdL);
    if((cntSmSigSmtdLDvert + cntSmBGSmtdLDvert) > 0)
	sgfSmtdLdVert = cntSmSigSmtdLDvert/sqrt(cntSmSigSmtdLDvert + cntSmBGSmtdLDvert);
    if((cntSmSigSmtdLDvertLm + cntSmBGSmtdLDvertLm) > 0)
	sgfSLcuts = cntSmSigSmtdLDvertLm/sqrt(cntSmSigSmtdLDvertLm + cntSmBGSmtdLDvertLm);
    if((cntSmSigSmtdS + cntSmBGSmtdS) > 0)
	sgfSLcutsMtdS = cntSmSigSmtdS/sqrt(cntSmSigSmtdS + cntSmBGSmtdS);

    cout << "significance: " << endl << "Lambda:" << " Lnc=" << sgfLnc << " LmtdL=" << sgfLmtdL << " LmtdLdVert=" << sgfLmtdLdVert <<  endl;
    cout << "Sigma: " << " Snc=" << sgfSnc << " SmtdL=" << sgfSmtdL << " SmtdLdVert=" << sgfSmtdLdVert <<
	" SLcuts=" << sgfSLcuts << " SLcutsMtdS=" << sgfSLcutsMtdS << endl;
    
   //reco effi
   double effLS, effLBG, effSS, effSBG;
   effLS = cntLmSigLmtdDvert/nentries*100;
//   effLBG = cntLmBGLmtdDlv/nentries*100;
   effSS = cntSmSigSmtdLDvertLm/nentries*100;
//   effSBG = cntSmBGSmtdLDlv/nentries*100;
   cout << "Effi: LMtdDvertLm = " << effLS << " SMtdLDvertSm = " << effSS  << endl;

     
   //edit histos
   hdEdx_Mdc -> GetXaxis() -> SetTitle("p*q [q*MeV/c]");
   hdEdx_Mdc -> GetYaxis() -> SetTitle("dE/dx [a.u.]");
   hdEdx_Mdc_acc -> GetXaxis() -> SetTitle("p*q [q*MeV/c]");
   hdEdx_Mdc_acc -> GetYaxis() -> SetTitle("dE/dx [a.u.]");

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
   hLDecVertX -> GetXaxis() -> SetTitle("#Lambda DecayVertX [mm]");
   hLDecVertY -> GetXaxis() -> SetTitle("#Lambda DecayVertY [mm]");
   hLDecVertZ -> GetXaxis() -> SetTitle("#Lambda DecayVertZ [mm]");
   hPrimVertR -> GetXaxis() -> SetTitle("PrimVertR [mm]");
   hLDecVertR -> GetXaxis() -> SetTitle("#Lambda DecayVertR [mm]");
   hdLVR -> GetXaxis() -> SetTitle("#Lambda DecVartR - PrimVertR [mm]");
   hdVert -> GetXaxis() -> SetTitle("dVert [mm]");
   hMinTrackDist -> GetXaxis() -> SetTitle("MTD_{p#pi^{-}} [mm]");
   hMinTrackDistLpi -> GetXaxis() -> SetTitle("MTD_{#Lambda#pi^{+}} [mm]");
   hMTD2D -> GetXaxis() -> SetTitle("MTD_{p#pi^{-}} [mm]");
   hMTD2D -> GetYaxis() -> SetTitle("MTD_{#Lambda#pi^{+}} [mm]");

   hLm -> GetYaxis() -> SetTitle("# [a.u]");
   hSm -> GetYaxis() -> SetTitle("# [a.u.]");
   hMinTrackDist -> GetYaxis() -> SetTitle("# [a.u.]");
   hMinTrackDistLpi -> GetYaxis() -> SetTitle("# [a.u.]");
   hPrimVertZ -> GetYaxis() -> SetTitle("# [a.u.]");
   hLDecVertX -> GetYaxis() -> SetTitle("# [a.u.]");
   hLDecVertY -> GetYaxis() -> SetTitle("# [a.u.]");
   hLDecVertZ -> GetYaxis() -> SetTitle("# [a.u.]");
   hPrimVertR -> GetYaxis() -> SetTitle("# [a.u.]");
   hLDecVertR -> GetYaxis() -> SetTitle("# [a.u.]");
   hdLVR -> GetYaxis() -> SetTitle("# [a.u.]");
   hdVert -> GetYaxis() -> SetTitle("# [a.u.]");

   
   hLm_mtdL -> GetXaxis() -> SetTitle("m_{p#pi^{-}} [MeV]");
   hSm_mtdL -> GetXaxis() -> SetTitle("m_{p#pi^{-}#pi^{+}} [MeV]");
   hLm_mtdL_dLV -> GetXaxis() -> SetTitle("m_{p#pi^{-}} [MeV]");
   hSm_mtdL_dLV -> GetXaxis() -> SetTitle("m_{p#pi^{-}#pi^{+}} [MeV]");
   hLm_mtdL_dvert -> GetXaxis() -> SetTitle("m_{p#pi^{-}} [MeV]");
   hSm_mtdL_dvert -> GetXaxis() -> SetTitle("m_{p#pi^{-}#pi^{+}} [MeV]");
   hLm_mtdL_dvert_Lm -> GetXaxis() -> SetTitle("m_{p#pi^{-}} [MeV]");
   hSm_mtdL_dvert_Lm -> GetXaxis() -> SetTitle("m_{p#pi^{-}#pi^{+}} [MeV]");
   hLm_mtdL_dLV_Lm -> GetXaxis() -> SetTitle("m_{p#pi^{-}} [MeV]");
   hSm_mtdL_dLV_Lm -> GetXaxis() -> SetTitle("m_{p#pi^{-}#pi^{+}} [MeV]");
   hLm_dLV -> GetXaxis() -> SetTitle("m_{p#pi^{-}} [MeV]");
   hSm_dLV -> GetXaxis() -> SetTitle("m_{p#pi^{-}#pi^{+}} [MeV]");
   hLm_dLV_Lm -> GetXaxis() -> SetTitle("m_{p#pi^{-}} [MeV]");
   hSm_dLV_Lm -> GetXaxis() -> SetTitle("m_{p#pi^{-}#pi^{+}} [MeV]");
   hSm_mtdS -> GetXaxis() -> SetTitle("m_{p#pi^{-}#pi^{+}} [MeV]");
   hSm_mtdS_pvz -> GetXaxis() -> SetTitle("m_{p#pi^{-}#pi^{+}} [MeV]");

   hLm_mtdL -> GetYaxis() -> SetTitle("# [a.u.]");
   hSm_mtdL -> GetYaxis() -> SetTitle("# [a.u.]");
   hLm_mtdL_dLV -> GetYaxis() -> SetTitle("# [a.u.]");
   hSm_mtdL_dLV -> GetYaxis() -> SetTitle("# [a.u.]");
   hLm_mtdL_dvert -> GetYaxis() -> SetTitle("# [a.u.]");
   hSm_mtdL_dvert -> GetYaxis() -> SetTitle("# [a.u.]");
   hLm_mtdL_dvert_Lm -> GetYaxis() -> SetTitle("# [a.u.]");
   hSm_mtdL_dvert_Lm -> GetYaxis() -> SetTitle("# [a.u.]");
   hLm_mtdL_dLV_Lm -> GetYaxis() -> SetTitle("# [a.u.]");
   hSm_mtdL_dLV_Lm -> GetYaxis() -> SetTitle("# [a.u.]");
   hLm_dLV -> GetYaxis() -> SetTitle("# [a.u.]");
   hSm_dLV -> GetYaxis() -> SetTitle("# [a.u.]");
   hLm_dLV_Lm -> GetYaxis() -> SetTitle("# [a.u.]");
   hSm_dLV_Lm -> GetYaxis() -> SetTitle("# [a.u.]");
   hSm_mtdS -> GetYaxis() -> SetTitle("# [a.u.]");
   hSm_mtdS_pvz -> GetYaxis() -> SetTitle("# [a.u.]");

   //   hLmRecBG -> GetXaxis() -> SetTitle("m_{p#pi^{-}} [MeV]");
   // histSig -> GetXaxis() -> SetTitle("m_{p#pi^{-}} [MeV]");
   // histBG -> GetXaxis() -> SetTitle("m_{p#pi^{-}} [MeV]");
   
   hMMABC -> GetXaxis() -> SetTitle("MM_{p#pi^{-}#pi^{+}} [MeV]");
   hMMABCpip -> GetXaxis() -> SetTitle("MM_{p#pi^{-}#pi^{+}} + M_{#pi^{+}} [MeV]");
   hMMABCpim -> GetXaxis() -> SetTitle("MM_{p#pi^{-}#pi^{+}} + M_{#pi^{-}} [MeV]");
   hinvMAC -> GetXaxis() -> SetTitle("invM_{p#pi^{+}} [MeV]");
   hinvMAB -> GetXaxis() -> SetTitle("invM_{p#pi^{-}} [MeV]");

   hDN -> GetXaxis() -> SetTitle("invM_{p#pi^{-}} [MeV]");
   hDN -> GetYaxis() -> SetTitle("MM_{p#pi^{-}#pi^{+}} + M_{#pi^{+}} [MeV]");
   hDpim -> GetXaxis() -> SetTitle("MM_{p#pi^{-}#pi^{+}} [MeV]");
   hDpim -> GetYaxis() -> SetTitle("invM_{p#pi^{-}} [MeV]");
   hDpip -> GetXaxis() -> SetTitle("MM_{p#pi^{-}#pi^{+}} [MeV]");
   hDpip -> GetYaxis() -> SetTitle("invM_{p#pi^{+}} [MeV]");

   hLm0 -> GetXaxis() -> SetTitle("m_{p#pi^{-}} [MeV]");
   hSm0 -> GetXaxis() -> SetTitle("m_{p#pi^{-}#pi^{+}} [MeV]");
   hMinTrackDist0 -> GetXaxis() -> SetTitle("MTD_{p#pi^{-}} [mm]");
   hMinTrackDistLpi0 -> GetXaxis() -> SetTitle("MTD_{#Lambda#pi^{+}} [mm]");
   hMTD2D0 -> GetXaxis() -> SetTitle("MTD_{p#pi^{-}} [mm]");
   hMTD2D0 -> GetYaxis() -> SetTitle("MTD_{#Lambda#pi^{+}} [mm]");
   hPrimVertZ0 -> GetXaxis() -> SetTitle("PrimVertZ [mm]");
   hLDecVertZ0 -> GetXaxis() -> SetTitle("#Lambda DecayVertZ [mm]");
   hPrimVertR0 -> GetXaxis() -> SetTitle("PrimVertR [mm]");
   hLDecVertR0 -> GetXaxis() -> SetTitle("#Lambda DecayVertR [mm]");
   hdLVR0 -> GetXaxis() -> SetTitle("#Lambda DecVartR - PrimVertR [mm]");
   hdVert0 -> GetXaxis() -> SetTitle("dVert [mm]");
   hMMABC_0 -> GetXaxis() -> SetTitle("MM_{p#pi^{-}#pi^{+}} [MeV]");
   hMMABC0 -> GetXaxis() -> SetTitle("MM_{p#pi^{-}#pi^{+}} [MeV]");
   hinvMAC0 -> GetXaxis() -> SetTitle("invM_{p#pi^{+}} [MeV]");
   hinvMAC2 -> GetXaxis() -> SetTitle("invM_{p#pi^{+}} [MeV]");
   hinvMAB0 -> GetXaxis() -> SetTitle("invM_{p#pi^{-}} [MeV]");
   hDpim0 -> GetXaxis() -> SetTitle("MM_{p#pi^{-}#pi^{+}} [MeV]");
   hDpip0 -> GetXaxis() -> SetTitle("MM_{p#pi^{-}#pi^{+}} [MeV]");

   hLm0 -> GetYaxis() -> SetTitle("# [a.u]");
   hSm0 -> GetYaxis() -> SetTitle("# [a.u.]");
   hMinTrackDist0 -> GetYaxis() -> SetTitle("# [a.u.]");
   hMinTrackDistLpi0 -> GetYaxis() -> SetTitle("# [a.u.]");
   hPrimVertZ0 -> GetYaxis() -> SetTitle("# [a.u.]");
   hLDecVertZ0 -> GetYaxis() -> SetTitle("# [a.u.]");
   hPrimVertR0 -> GetYaxis() -> SetTitle("# [a.u.]");
   hLDecVertR0 -> GetYaxis() -> SetTitle("# [a.u.]");
   hdLVR0 -> GetYaxis() -> SetTitle("# [a.u.]");
   hdVert0 -> GetYaxis() -> SetTitle("# [a.u.]");
   hMMABC_0 -> GetYaxis() -> SetTitle("# [a.u.]");
   hMMABC0 -> GetYaxis() -> SetTitle("# [a.u.]");
   hinvMAC0 -> GetYaxis() -> SetTitle("# [a.u.]");
   hinvMAC2 -> GetYaxis() -> SetTitle("# [a.u.]");
   hinvMAB0 -> GetYaxis() -> SetTitle("# [a.u.]");
   hDpim0 -> GetYaxis() -> SetTitle("invM_{p#pi^{-}} [MeV]");
   hDpip0 -> GetYaxis() -> SetTitle("invM_{p#pi^{+}} [MeV]");

   TLine *limm1 = new TLine(cmm1, 0, cmm1, dmaxy);
   TLine *limm2 = new TLine(0, cmm2, dmaxx, cmm2);
   limm1 -> SetLineWidth(2);
   limm2 -> SetLineWidth(2);
   cDpip -> cd();
   gPad -> SetLogz();
   hDpip0 -> Draw("colz");
   limm1 -> Draw("same");
   limm2 -> Draw("same");
   
   //vert_sim - vert_rec
/*   hdvertlX -> GetXaxis() -> SetTitle("d(vertL_sim - vertL_rec)x [mm]");
   hdvertlY -> GetXaxis() -> SetTitle("d(vertL_sim - vertL_rec)y [mm]");
   hdvertlZ -> GetXaxis() -> SetTitle("d(vertL_sim - vertL_rec)z [mm]");
   hdvertlR -> GetXaxis() -> SetTitle("d(vertL_sim - vertL_rec)r [mm]");
   hdvertlX -> GetYaxis() -> SetTitle("counts");
   hdvertlY -> GetYaxis() -> SetTitle("counts");
   hdvertlZ -> GetYaxis() -> SetTitle("counts");
   hdvertlR -> GetYaxis() -> SetTitle("counts");
   TCanvas *cdvertl = new TCanvas("cdvertl", "cdvertl");
   cdvertl -> Divide(2,2);
   cdvertl->cd(1);
   hdvertlX -> Draw();
   cdvertl->cd(2);
   hdvertlY -> Draw();
   cdvertl->cd(3);
   hdvertlZ -> Draw();
   cdvertl->cd(4);
   hdvertlR -> Draw();

   //vert_rec
   hvx_lambdaX -> GetXaxis() -> SetTitle("(vertL_rec)x [mm]");
   hvx_lambdaY -> GetXaxis() -> SetTitle("(vertL_rec)y [mm]");
   hvx_lambdaZ -> GetXaxis() -> SetTitle("(vertL_rec)z [mm]");
   hvx_lambdaX -> GetYaxis() -> SetTitle("counts");
   hvx_lambdaY -> GetYaxis() -> SetTitle("counts");
   hvx_lambdaZ -> GetYaxis() -> SetTitle("counts");
   TCanvas *cvx_lambda = new TCanvas("cvx_lambda", "cvx_lambda");
   cvx_lambda -> Divide(2,2);
   cvx_lambda->cd(1);
   hvx_lambdaX -> Draw();
   cvx_lambda->cd(2);
   hvx_lambdaY -> Draw();
   cvx_lambda->cd(3);
   hvx_lambdaZ -> Draw();

   //vert_sim
   hgeantxvertexA -> GetXaxis() -> SetTitle("(vertL_sim)x [mm]");
   hgeantyvertexA -> GetXaxis() -> SetTitle("(vertL_sim)y [mm]");
   hgeantzvertexA -> GetXaxis() -> SetTitle("(vertL_sim)z [mm]");
   hgeantxvertexA -> GetYaxis() -> SetTitle("counts");
   hgeantyvertexA -> GetYaxis() -> SetTitle("counts");
   hgeantzvertexA -> GetYaxis() -> SetTitle("counts");
   TCanvas *cgeantvertexA = new TCanvas("cgeantvertexA", "cgeantvertexA");
   cgeantvertexA -> Divide(2,2);
   cgeantvertexA->cd(1);
   hgeantxvertexA -> Draw();
   cgeantvertexA->cd(2);
   hgeantyvertexA -> Draw();
   cgeantvertexA->cd(3);
   hgeantzvertexA -> Draw();

   //vert_sim_reco
   hvertex_recoGeant_x -> GetXaxis() -> SetTitle("(vertL_sim_reco)x [mm]");
   hvertex_recoGeant_y -> GetXaxis() -> SetTitle("(vertL_sim_reco)y [mm]");
   hvertex_recoGeant_z -> GetXaxis() -> SetTitle("(vertL_sim_reco)z [mm]");
   hvertex_recoGeant_x -> GetYaxis() -> SetTitle("counts");
   hvertex_recoGeant_y -> GetYaxis() -> SetTitle("counts");
   hvertex_recoGeant_z -> GetYaxis() -> SetTitle("counts");
   TCanvas *cvertex_recoGeant = new TCanvas("cvertex_recoGeant", "cvertex_recoGeant");
   cvertex_recoGeant -> Divide(2,2);
   cvertex_recoGeant -> cd(1);
   hvertex_recoGeant_x -> Draw();
   cvertex_recoGeant -> cd(2);
   hvertex_recoGeant_y -> Draw();
   cvertex_recoGeant -> cd(3);
   hvertex_recoGeant_z -> Draw();
   
   //proton momentum from simulation
   hMomAGeax -> GetXaxis() -> SetTitle("(momPsim)x ");
   hMomAGeay -> GetXaxis() -> SetTitle("(momPsim)y ");
   hMomAGeaz -> GetXaxis() -> SetTitle("(momPsim)z ");
   hMomAGeax -> GetYaxis() -> SetTitle("counts");
   hMomAGeay -> GetYaxis() -> SetTitle("counts");
   hMomAGeaz -> GetYaxis() -> SetTitle("counts");
   TCanvas *cMomAGea = new TCanvas("cMomAGea", "cMomAGea");
   cMomAGea -> Divide(2,2);
   cMomAGea -> cd(1);
   hMomAGeax -> Draw();
   cMomAGea -> cd(2);
   hMomAGeay -> Draw();
   cMomAGea -> cd(3);
   hMomAGeaz -> Draw();
   
   //proton momentum reconstructed
   hMomArecox -> GetXaxis() -> SetTitle("(momPreco)x ");
   hMomArecoy -> GetXaxis() -> SetTitle("(momPreco)y ");
   hMomArecoz -> GetXaxis() -> SetTitle("(momPreco)z ");
   hMomArecox -> GetYaxis() -> SetTitle("counts");
   hMomArecoy -> GetYaxis() -> SetTitle("counts");
   hMomArecoz -> GetYaxis() -> SetTitle("counts");
   TCanvas *cMomAreco = new TCanvas("cMomAreco", "cMomAreco");
   cMomAreco -> Divide(2,2);
   cMomAreco -> cd(1);
   hMomArecox -> Draw();
   cMomAreco -> cd(2);
   hMomArecoy -> Draw();
   cMomAreco -> cd(3);
   hMomArecoz -> Draw();
   
   //proton momentum resolution
   hMomAResx -> GetXaxis() -> SetTitle("(momPreco - momPsim)x ");
   hMomAResy -> GetXaxis() -> SetTitle("(momPreco - momPsim)y ");
   hMomAResz -> GetXaxis() -> SetTitle("(momPreco - momPsim)z ");
   hMomAResx -> GetYaxis() -> SetTitle("counts");
   hMomAResy -> GetYaxis() -> SetTitle("counts");
   hMomAResz -> GetYaxis() -> SetTitle("counts");
   TCanvas *cMomARes = new TCanvas("cMomARes", "cMomARes");
   cMomARes -> Divide(2,2);
   cMomARes -> cd(1);
   hMomAResx -> Draw();
   cMomARes -> cd(2);
   hMomAResy -> Draw();
   cMomARes -> cd(3);
   hMomAResz -> Draw();

   //L momentum from simulation
   hMomGeax -> GetXaxis() -> SetTitle("(momLsim)x ");
   hMomGeay -> GetXaxis() -> SetTitle("(momLsim)y ");
   hMomGeaz -> GetXaxis() -> SetTitle("(momLsim)z ");
   hMomGeax -> GetYaxis() -> SetTitle("counts");
   hMomGeay -> GetYaxis() -> SetTitle("counts");
   hMomGeaz -> GetYaxis() -> SetTitle("counts");
   TCanvas *cMomGea = new TCanvas("cMomGea", "cMomGea");
   cMomGea -> Divide(2,2);
   cMomGea -> cd(1);
   hMomGeax -> Draw();
   cMomGea -> cd(2);
   hMomGeay -> Draw();
   cMomGea -> cd(3);
   hMomGeaz -> Draw();
   
   //L momentum reconstructed
   hMomxLreco -> GetXaxis() -> SetTitle("(momLreco)x ");
   hMomyLreco -> GetXaxis() -> SetTitle("(momLreco)y ");
   hMomzLreco -> GetXaxis() -> SetTitle("(momLreco)z ");
   hMomxLreco -> GetYaxis() -> SetTitle("counts");
   hMomyLreco -> GetYaxis() -> SetTitle("counts");
   hMomzLreco -> GetYaxis() -> SetTitle("counts");
   TCanvas *cMomLreco = new TCanvas("cMomLreco", "cMomLreco");
   cMomLreco -> Divide(2,2);
   cMomLreco -> cd(1);
   hMomxLreco -> Draw();
   cMomLreco -> cd(2);
   hMomyLreco -> Draw();
   cMomLreco -> cd(3);
   hMomzLreco -> Draw();
   
   //L momentum resolution
   hMomResx -> GetXaxis() -> SetTitle("(momLreco - momLsim)x ");
   hMomResy -> GetXaxis() -> SetTitle("(momLreco - momLsim)y ");
   hMomResz -> GetXaxis() -> SetTitle("(momLreco - momLsim)z ");
   hMomResx -> GetYaxis() -> SetTitle("counts");
   hMomResy -> GetYaxis() -> SetTitle("counts");
   hMomResz -> GetYaxis() -> SetTitle("counts");
   TCanvas *cMomRes = new TCanvas("cMomRes", "cMomRes");
   cMomRes -> Divide(2,2);
   cMomRes -> cd(1);
   hMomResx -> Draw();
   cMomRes -> cd(2);
   hMomResy -> Draw();
   cMomRes -> cd(3);
   hMomResz -> Draw();

   hLmass2 -> GetXaxis() -> SetTitle("L mass [MeV]");
   hLmass2 -> GetYaxis() -> SetTitle("counts");
*/   
   /* hrA -> GetXaxis() -> SetTitle("r_#pi^{-} [mm]");
   hrB -> GetXaxis() -> SetTitle("r_p [mm]");
   hrA -> GetYaxis() -> SetTitle("counts");
   hrB -> GetYaxis() -> SetTitle("counts");

   hrealLambdaM -> GetXaxis() -> SetTitle("m(#Lambda) [meV]");
   hrealLambdaM -> GetYaxis() -> SetTitle("counts");
   */
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
   cLsbg1 -> cd(1);
   nice_canv1(gPad);
   hLm0 -> Draw();
   fhistLnc0 -> Draw("same");
   fbgLnc0 -> Draw("same p");
   fsigLnc0 -> Draw("same p");
//   histSigLnc0 -> Draw("same");
   cLsbg1 -> cd(2);
   nice_canv1(gPad);
   hLm -> Draw();
   fhistLnc -> Draw("same");
   fbgLnc -> Draw("same p");
   fsigLnc -> Draw("same p");
//   histSigLnc -> Draw("same");
   cLsbg2 -> cd(1);
   nice_canv1(gPad);
   hLm_mtdL -> Draw();
   fhistLmtd -> Draw("same");
   fbgLmtd -> Draw("same p");
   fsigLmtd -> Draw("same p");
//   histSigLmtd -> Draw("same");
   cLsbg2 -> cd(2);
   nice_canv1(gPad);
   hLm_mtdL_dvert -> Draw();
   fhistLmtdDvert -> Draw("same");
   fbgLmtdDvert -> Draw("same p");
   fsigLmtdDvert -> Draw("same p");
//   histSigLmtdDlv -> Draw("same");
   //S
   cSsbg1 -> cd(1);
   nice_canv1(gPad);
   hSm0 -> Draw();
   //fbgSnc0 -> Draw("same p");
   //histSigSnc0 -> Draw("same");
   cSsbg1 -> cd(2);
   nice_canv1(gPad);
   hSm -> Draw();
   //fbgSnc -> Draw("same p");
   //histSigSnc -> Draw("same");
   cSsbg2 -> cd(1);
   nice_canv1(gPad);
   hSm_mtdL -> Draw();
   //fsigSmtdL -> Draw("same p");
   //histSigSmtdL -> Draw("same");
   cSsbg2 -> cd(2);
   nice_canv1(gPad);
   hSm_mtdL_dvert -> Draw();
   //fsigSmtdLDlv -> Draw("same p");
   //histSigSmtdLDlv -> Draw("same");
   cSsbg2 -> cd(3);
   nice_canv1(gPad);
   hSm_mtdL_dvert_Lm -> Draw();
   fhistSmtdLDvertLm -> Draw("same");
   fsigSmtdLDvertLm -> Draw("same p");
   fbgSmtdLDvertLm -> Draw("same p");
   ptBW -> Draw("same");
   
   
   //writing histos
   TFile *fout = TFile::Open("./outputs/anaSigmaOut_exp_all0919_3.root", "RECREATE");
//   TFile *fout = TFile::Open("./outputs/test.root", "RECREATE");
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
//<<<<
/*
   cMtdScan -> Write();
   cMtdScanSigma -> Write();
   gmtdscan -> Write();
   //cdLVScan -> Write();
   gdVertscan -> Write();
   //gdLVmtdL -> Write();
//   gmtdLpiscan -> Write();
   gdVertmtdL -> Write();
   hdVertmtdL -> Write();
   for(int i = 0; i < 15; i++){
//       cmtdLdLV[i] -> Write();
       cmtdLdvert[i] -> Write();
       cmtdLdvertSigma[i] -> Write();       
   }
//   cMtdLpiScan -> Write();
//<<<<
   */

//
   hPrimVertZ -> Write();
   hLDecVertX -> Write();
   hLDecVertY -> Write();
   hLDecVertZ -> Write();
   hPrimVertR -> Write();
   hLDecVertR -> Write();
   hdLVR -> Write();
   hdVert -> Write();
   hMinTrackDist -> Write();
   hMinTrackDistLpi -> Write();
   hMTD2D -> Write();
   
   hLm_mtdL -> Write();
   hSm_mtdL -> Write();
   hLm_mtdL_dvert -> Write();
   hSm_mtdL_dvert -> Write();
   hLm_mtdL_dvert_Lm -> Write();
   hSm_mtdL_dvert_Lm -> Write();
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

   cLsbg1 -> Write();
   cLsbg2 -> Write();
   cSsbg1 -> Write();
   cSsbg2 -> Write();
   
   hMMABC -> Write();
   hMMABCpip -> Write();
   hMMABCpim -> Write();
   hinvMAC -> Write();
   hinvMAB -> Write();
   hMMAB -> Write();
   hMMBC -> Write();
   hDN -> Write();
   hDpim -> Write();

   hDpip -> Write();
   hLm0 -> Write();
   hSm0 -> Write();
   hMinTrackDist0 -> Write();
   hMinTrackDistLpi0 -> Write();
   hMTD2D0 -> Write();
   hPrimVertR0 -> Write();
   hLDecVertR0 -> Write();
   hdLVR0 -> Write();
   hdVert0 -> Write();
   hPrimVertZ0 -> Write();
   hLDecVertZ0 -> Write();
   hMMABC_0 -> Write();
   hMMABC0 -> Write();
   hinvMAC0 -> Write();
   hinvMAC2 -> Write();
   hinvMAB0 -> Write();
   hDpim0 -> Write();
   hDpip0 -> Write();
   cDpip -> Write();

   cSigmaFit ->Write();

/*   hdvertlX -> Write();
   hdvertlY -> Write();
   hdvertlZ -> Write();
   hdvertlR -> Write();
   cdvertl -> Write();
     
   hvx_lambdaX -> Write();
   hvx_lambdaY -> Write();
   hvx_lambdaZ -> Write();
   cvx_lambda -> Write();
     
   hgeantxvertexA -> Write();
   hgeantyvertexA -> Write();
   hgeantzvertexA -> Write();
   cgeantvertexA -> Write();
     
   hvertex_recoGeant_x -> Write();
   hvertex_recoGeant_y -> Write();
   hvertex_recoGeant_z -> Write();
   cvertex_recoGeant -> Write();

   hMomAGeax -> Write();
   hMomAGeay -> Write();
   hMomAGeaz -> Write();
   cMomAGea -> Write();
   hMomArecox -> Write();
   hMomArecoy -> Write();
   hMomArecoz -> Write();
   cMomAreco -> Write();
   hMomAResx -> Write();
   hMomAResy -> Write();
   hMomAResz -> Write();
   cMomARes -> Write();

   hMomGeax -> Write();
   hMomGeay -> Write();
   hMomGeaz -> Write();
   cMomGea -> Write();
   hMomxLreco -> Write();
   hMomyLreco -> Write();
   hMomzLreco -> Write();
   cMomLreco -> Write();
   hMomResx -> Write();
   hMomResy -> Write();
   hMomResz -> Write();
   cMomRes -> Write();

   hLmass2 -> Write();
*/   
   /* hrA -> Write();
   hrB -> Write();
   
   hrealLambdaM -> Write();
   */   
   heventNo -> Write();
   
   fout -> Close();
}
