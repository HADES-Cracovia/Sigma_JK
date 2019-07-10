#include <TH1F.h>

double peakBG(TH1F* h)
{
//=========Macro generated from canvas: Canvas_1/Canvas_1
//=========  (Wed Jul 10 11:07:50 2019) by ROOT version 6.14/06
    /* TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",258,116,1658,838);
   Canvas_1->ToggleEventStatus();
   Canvas_1->Range(1100,-2444.663,2100,22001.96);
   Canvas_1->SetFillColor(0);
   Canvas_1->SetBorderMode(0);
   Canvas_1->SetBorderSize(2);
   Canvas_1->SetFrameBorderMode(0);
   Canvas_1->SetFrameBorderMode(0);
    */
    TH1F *hSm_mtdL_dvert_Lm__1 = (TH1F*)h -> Clone();
/*   hSm_mtdL_dvert_Lm__1->SetBinContent(12,28);
   hSm_mtdL_dvert_Lm__1->SetBinContent(13,266);
   hSm_mtdL_dvert_Lm__1->SetBinContent(14,744);
   hSm_mtdL_dvert_Lm__1->SetBinContent(15,1431);
   hSm_mtdL_dvert_Lm__1->SetBinContent(16,2272);
   hSm_mtdL_dvert_Lm__1->SetBinContent(17,3309);
   hSm_mtdL_dvert_Lm__1->SetBinContent(18,4410);
   hSm_mtdL_dvert_Lm__1->SetBinContent(19,5455);
   hSm_mtdL_dvert_Lm__1->SetBinContent(20,6129);
   hSm_mtdL_dvert_Lm__1->SetBinContent(21,7147);
   hSm_mtdL_dvert_Lm__1->SetBinContent(22,7985);
   hSm_mtdL_dvert_Lm__1->SetBinContent(23,8712);
   hSm_mtdL_dvert_Lm__1->SetBinContent(24,9328);
   hSm_mtdL_dvert_Lm__1->SetBinContent(25,9949);
   hSm_mtdL_dvert_Lm__1->SetBinContent(26,10468);
   hSm_mtdL_dvert_Lm__1->SetBinContent(27,11052);
   hSm_mtdL_dvert_Lm__1->SetBinContent(28,11734);
   hSm_mtdL_dvert_Lm__1->SetBinContent(29,12272);
   hSm_mtdL_dvert_Lm__1->SetBinContent(30,12680);
   hSm_mtdL_dvert_Lm__1->SetBinContent(31,13098);
   hSm_mtdL_dvert_Lm__1->SetBinContent(32,13714);
   hSm_mtdL_dvert_Lm__1->SetBinContent(33,14114);
   hSm_mtdL_dvert_Lm__1->SetBinContent(34,14768);
   hSm_mtdL_dvert_Lm__1->SetBinContent(35,15061);
   hSm_mtdL_dvert_Lm__1->SetBinContent(36,15624);
   hSm_mtdL_dvert_Lm__1->SetBinContent(37,15957);
   hSm_mtdL_dvert_Lm__1->SetBinContent(38,16575);
   hSm_mtdL_dvert_Lm__1->SetBinContent(39,17035);
   hSm_mtdL_dvert_Lm__1->SetBinContent(40,17150);
   hSm_mtdL_dvert_Lm__1->SetBinContent(41,17496);
   hSm_mtdL_dvert_Lm__1->SetBinContent(42,17902);
   hSm_mtdL_dvert_Lm__1->SetBinContent(43,18185);
   hSm_mtdL_dvert_Lm__1->SetBinContent(44,18199);
   hSm_mtdL_dvert_Lm__1->SetBinContent(45,18626);
   hSm_mtdL_dvert_Lm__1->SetBinContent(46,18414);
   hSm_mtdL_dvert_Lm__1->SetBinContent(47,17958);
   hSm_mtdL_dvert_Lm__1->SetBinContent(48,17798);
   hSm_mtdL_dvert_Lm__1->SetBinContent(49,17403);
   hSm_mtdL_dvert_Lm__1->SetBinContent(50,17129);
   hSm_mtdL_dvert_Lm__1->SetBinContent(51,16664);
   hSm_mtdL_dvert_Lm__1->SetBinContent(52,16157);
   hSm_mtdL_dvert_Lm__1->SetBinContent(53,16111);
   hSm_mtdL_dvert_Lm__1->SetBinContent(54,15825);
   hSm_mtdL_dvert_Lm__1->SetBinContent(55,15455);
   hSm_mtdL_dvert_Lm__1->SetBinContent(56,15392);
   hSm_mtdL_dvert_Lm__1->SetBinContent(57,15566);
   hSm_mtdL_dvert_Lm__1->SetBinContent(58,15583);
   hSm_mtdL_dvert_Lm__1->SetBinContent(59,15817);
   hSm_mtdL_dvert_Lm__1->SetBinContent(60,15962);
   hSm_mtdL_dvert_Lm__1->SetBinContent(61,16053);
   hSm_mtdL_dvert_Lm__1->SetBinContent(62,16439);
   hSm_mtdL_dvert_Lm__1->SetBinContent(63,15794);
   hSm_mtdL_dvert_Lm__1->SetBinContent(64,15817);
   hSm_mtdL_dvert_Lm__1->SetBinContent(65,15791);
   hSm_mtdL_dvert_Lm__1->SetBinContent(66,15476);
   hSm_mtdL_dvert_Lm__1->SetBinContent(67,14803);
   hSm_mtdL_dvert_Lm__1->SetBinContent(68,14557);
   hSm_mtdL_dvert_Lm__1->SetBinContent(69,13888);
   hSm_mtdL_dvert_Lm__1->SetBinContent(70,13449);
   hSm_mtdL_dvert_Lm__1->SetBinContent(71,12519);
   hSm_mtdL_dvert_Lm__1->SetBinContent(72,11613);
   hSm_mtdL_dvert_Lm__1->SetBinContent(73,10614);
   hSm_mtdL_dvert_Lm__1->SetBinContent(74,9654);
   hSm_mtdL_dvert_Lm__1->SetBinContent(75,8772);
   hSm_mtdL_dvert_Lm__1->SetBinContent(76,7922);
   hSm_mtdL_dvert_Lm__1->SetBinContent(77,7595);
   hSm_mtdL_dvert_Lm__1->SetBinContent(78,7223);
   hSm_mtdL_dvert_Lm__1->SetBinContent(79,6889);
   hSm_mtdL_dvert_Lm__1->SetBinContent(80,6637);
   hSm_mtdL_dvert_Lm__1->SetBinContent(81,6481);
   hSm_mtdL_dvert_Lm__1->SetBinContent(82,6215);
   hSm_mtdL_dvert_Lm__1->SetBinContent(83,6119);
   hSm_mtdL_dvert_Lm__1->SetBinContent(84,5925);
   hSm_mtdL_dvert_Lm__1->SetBinContent(85,5781);
   hSm_mtdL_dvert_Lm__1->SetBinContent(86,5617);
   hSm_mtdL_dvert_Lm__1->SetBinContent(87,5526);
   hSm_mtdL_dvert_Lm__1->SetBinContent(88,5272);
   hSm_mtdL_dvert_Lm__1->SetBinContent(89,5000);
   hSm_mtdL_dvert_Lm__1->SetBinContent(90,4961);
   hSm_mtdL_dvert_Lm__1->SetBinContent(91,4701);
   hSm_mtdL_dvert_Lm__1->SetBinContent(92,4441);
   hSm_mtdL_dvert_Lm__1->SetBinContent(93,4479);
   hSm_mtdL_dvert_Lm__1->SetBinContent(94,4330);
   hSm_mtdL_dvert_Lm__1->SetBinContent(95,4163);
   hSm_mtdL_dvert_Lm__1->SetBinContent(96,3961);
   hSm_mtdL_dvert_Lm__1->SetBinContent(97,3865);
   hSm_mtdL_dvert_Lm__1->SetBinContent(98,3687);
   hSm_mtdL_dvert_Lm__1->SetBinContent(99,3541);
   hSm_mtdL_dvert_Lm__1->SetBinContent(100,3567);
   hSm_mtdL_dvert_Lm__1->SetBinContent(101,3354);
   hSm_mtdL_dvert_Lm__1->SetBinContent(102,3222);
   hSm_mtdL_dvert_Lm__1->SetBinContent(103,2992);
   hSm_mtdL_dvert_Lm__1->SetBinContent(104,2861);
   hSm_mtdL_dvert_Lm__1->SetBinContent(105,2843);
   hSm_mtdL_dvert_Lm__1->SetBinContent(106,2688);
   hSm_mtdL_dvert_Lm__1->SetBinContent(107,2504);
   hSm_mtdL_dvert_Lm__1->SetBinContent(108,2392);
   hSm_mtdL_dvert_Lm__1->SetBinContent(109,2275);
   hSm_mtdL_dvert_Lm__1->SetBinContent(110,2223);
   hSm_mtdL_dvert_Lm__1->SetBinContent(111,2069);
   hSm_mtdL_dvert_Lm__1->SetBinContent(112,2021);
   hSm_mtdL_dvert_Lm__1->SetBinContent(113,1954);
   hSm_mtdL_dvert_Lm__1->SetBinContent(114,1821);
   hSm_mtdL_dvert_Lm__1->SetBinContent(115,1725);
   hSm_mtdL_dvert_Lm__1->SetBinContent(116,1493);
   hSm_mtdL_dvert_Lm__1->SetBinContent(117,1465);
   hSm_mtdL_dvert_Lm__1->SetBinContent(118,1359);
   hSm_mtdL_dvert_Lm__1->SetBinContent(119,1337);
   hSm_mtdL_dvert_Lm__1->SetBinContent(120,1284);
   hSm_mtdL_dvert_Lm__1->SetBinContent(121,1114);
   hSm_mtdL_dvert_Lm__1->SetBinContent(122,1034);
   hSm_mtdL_dvert_Lm__1->SetBinContent(123,1030);
   hSm_mtdL_dvert_Lm__1->SetBinContent(124,920);
   hSm_mtdL_dvert_Lm__1->SetBinContent(125,842);
   hSm_mtdL_dvert_Lm__1->SetBinContent(126,786);
   hSm_mtdL_dvert_Lm__1->SetBinContent(127,733);
   hSm_mtdL_dvert_Lm__1->SetBinContent(128,700);
   hSm_mtdL_dvert_Lm__1->SetBinContent(129,637);
   hSm_mtdL_dvert_Lm__1->SetBinContent(130,585);
   hSm_mtdL_dvert_Lm__1->SetBinContent(131,510);
   hSm_mtdL_dvert_Lm__1->SetBinContent(132,458);
   hSm_mtdL_dvert_Lm__1->SetBinContent(133,416);
   hSm_mtdL_dvert_Lm__1->SetBinContent(134,370);
   hSm_mtdL_dvert_Lm__1->SetBinContent(135,322);
   hSm_mtdL_dvert_Lm__1->SetBinContent(136,346);
   hSm_mtdL_dvert_Lm__1->SetBinContent(137,255);
   hSm_mtdL_dvert_Lm__1->SetBinContent(138,269);
   hSm_mtdL_dvert_Lm__1->SetBinContent(139,217);
   hSm_mtdL_dvert_Lm__1->SetBinContent(140,188);
   hSm_mtdL_dvert_Lm__1->SetBinContent(141,187);
   hSm_mtdL_dvert_Lm__1->SetBinContent(142,146);
   hSm_mtdL_dvert_Lm__1->SetBinContent(143,145);
   hSm_mtdL_dvert_Lm__1->SetBinContent(144,103);
   hSm_mtdL_dvert_Lm__1->SetBinContent(145,99);
   hSm_mtdL_dvert_Lm__1->SetBinContent(146,85);
   hSm_mtdL_dvert_Lm__1->SetBinContent(147,58);
   hSm_mtdL_dvert_Lm__1->SetBinContent(148,50);
   hSm_mtdL_dvert_Lm__1->SetBinContent(149,43);
   hSm_mtdL_dvert_Lm__1->SetBinContent(150,49);
   hSm_mtdL_dvert_Lm__1->SetBinContent(151,32);
   hSm_mtdL_dvert_Lm__1->SetBinContent(152,13);
   hSm_mtdL_dvert_Lm__1->SetBinContent(153,12);
   hSm_mtdL_dvert_Lm__1->SetBinContent(154,6);
   hSm_mtdL_dvert_Lm__1->SetBinContent(155,8);
   hSm_mtdL_dvert_Lm__1->SetBinContent(156,3);
   hSm_mtdL_dvert_Lm__1->SetEntries(1011899);
*/
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
//       cout << " l x y: " << l << " " <<  x[l] << " " << y[l] << endl;
       l++;
     }
     j++;
   }

   TGraph *gr = new TGraph(n, x, y);
   gr -> Draw("ap*");

   int fitS1 = 1255;
   int fitS2 = 1600;
   int fitSsig1 = 1325;
   int fitSsig2 = 1445;

/*   TF1 * f1 = new TF1("f1", "[0]*[1]*[1]/( (x*x-[2]*[2])*(x*x-[2]*[2]) + (x*x*x*x*[1]*[1]/([2]*[2])) ) + pol6(3)", fitS1, fitS2);
   f1 -> SetParameters(
       1.40794e+09, 36, 1.38415e+03,
       -9.08964e+05,2.35372e+03,-2.33096e+00,
       1.08507e-03,-2.28651e-07,1.57016e-11
       );
   f1 -> SetParLimits(0,0,10000000000);
   f1 -> SetParLimits(1,30,40);
   f1 -> SetParLimits(2,1360,1400);

   gr -> Fit("f1", "0", "", fitS1, fitS2);
   gr -> Fit("f1", "0", "", fitS1, fitS2);

   f1 -> SetLineColor(kGreen);
*/ 
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
   
//   fhistSmtdLDvertLm -> Draw("same");
//   fsigSmtdLDvertLm -> Draw("same");
//   fbgSmtdLDvertLm -> Draw("same");
//   f1 -> Draw("same");
//   hSm_mtdL_dvert_Lm__1->Draw("same");

   return parSmtdLDvertLm[];
}
