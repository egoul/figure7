#define analysis_cxx
#include "analysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <cmath>
#include <math.h>
#include "makeMultiPanelCanvas.c"

// this code attempts to recreate the distribution of jetpt/gammapt analagous to figure 7 in the paper http://cds.cern.ch/record/2217884/files/HIN-16-002-pas.pdf
float allpuritypbpb[] = {0.725758, 0.720249, 0.753094, 0.703853, 0.730487, 0.756007, 0.741809, 0.737945, 0.725995, 0.786819, 0.704426, 0.743147, 0.775786, 0.827101, 0.715906, 0.710863, 0.739059, 0.687001, 0.719965, 0.743805, 0.720358, 0.724734, 0.719121, 0.758637, 0.699659, 0.73539, 0.76691, 0.731031, 0.730948, 0.717619, 0.786722, 0.695959, 0.733827, 0.771255, 0.836782, 0.749695, 0.739283, 0.788322, 0.711966, 0.763817, 0.786493, 0.816522, 0.733207, 0.719114, 0.785248, 0.705784, 0.734519, 0.769831, 0.84767, 0.772594, 0.766181, 0.802428, 0.720997, 0.80606, 0.822834, 0.774406};
float allpuritypp[] = {0.823368, 0.823368, 0.823368, 0.823368, 0.823368, 0.823368, 0.823368, 0.846154, 0.846154, 0.846154, 0.846154, 0.846154, 0.846154, 0.846154, 0.820975, 0.820975, 0.820975, 0.820975, 0.820975, 0.820975, 0.820975, 0.830048, 0.830048, 0.830048, 0.830048, 0.830048, 0.830048, 0.830048, 0.846293, 0.846293, 0.846293, 0.846293, 0.846293, 0.846293, 0.846293, 0.859037, 0.859037, 0.859037, 0.859037, 0.859037, 0.859037, 0.859037, 0.863744, 0.863744, 0.863744, 0.863744, 0.863744, 0.863744, 0.863744, 0.857244, 0.857244, 0.857244, 0.857244, 0.857244, 0.857244, 0.857244};
float getpurity(float phoetmin, float hibinmin, bool ispp)
{
    int row = -1;
    int col = -1;
    if(phoetmin==40)  row = 0;
    if(phoetmin==60)  row = 1;
    if(phoetmin==100) row = 7;
    if(hibinmin==0)   col = 3;
    if(hibinmin==20)  col = 4;
    if(hibinmin==60)  col = 5;
    if(hibinmin==100) col = 6;
    if(row>-1 && col > -1 && ispp) return allpuritypp[row*7+col];
    if(row>-1 && col > -1 && !ispp) return allpuritypbpb[row*7+col];
    return 1; //no purity applied
}

void analysis::Loop()
{
//   In a ROOT session, you can do:
//      root> .L analysis.C
//      root> analysis t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
    //TCanvas * call = new TCanvas("call","",1600,500);
    //makeMultiPanelCanvas(call,ncentbins+1,1,0.02,0.0,-6,0.2,0.04);
    
    TFile * outputfile = new TFile("meantrkPt.root","recreate");
    TH1D * h1 = new TH1D("h1test","title;xaxis title;yaxis title",20,0,3.4);
    TH1D * h1m = new TH1D("h1mixed","title;xaxis title;yaxis title",20,0,3.4);
    TH1D * h2 = new TH1D("h2test","title;xaxis title;yaxis title",20,0,3);
    TH1D * h2m = new TH1D("h2mixed","title;xaxis title;yaxis title",20,0,3);
    TH1D * h3 = new TH1D("h3test","title;xaxis title;yaxis title",16,0,2);
    TH1D * h3m = new TH1D("h3mixed","title;xaxis title;yaxis title",16,0,2);
    TH1D * h4 = new TH1D("h4test","title;xaxis title;yaxis title",20,0,3.4);
    TH1D * h4m = new TH1D("h4mixed","title;xaxis title;yaxis title",20,0,3.4);
    
    TH1D * h1sig = new TH1D("h1sig","title;xaxis title;yaxis title",20,0,3.4);
    TH1D * h1bg = new TH1D("h1bg","title;xaxis title;yaxis title",20,0,3.4);
    TH1D * h2sig = new TH1D("h2sig","title;xaxis title;yaxis title",20,0,3.4);
    TH1D * h2bg = new TH1D("h2bg","title;xaxis title;yaxis title",20,0,3.4);
    TH1D * h3sig = new TH1D("h3sig","title;xaxis title;yaxis title",20,0,3.4);
    TH1D * h3bg = new TH1D("h3bg","title;xaxis title;yaxis title",20,0,3.4);
    TH1D * h4sig = new TH1D("h4sig","title;xaxis title;yaxis title",20,0,3.4);
    TH1D * h4bg = new TH1D("h4bg","title;xaxis title;yaxis title",20,0,3.4);
    
    TH1D * h1f = new TH1D("h1final","title;xaxis title;yaxis title",20,0,3.4);
    TH1D * h2f = new TH1D("h2final","title;xaxis title;yaxis title",20,0,3.4);
    TH1D * h3f = new TH1D("h3final","title;xaxis title;yaxis title",20,0,3.4);
    TH1D * h4f = new TH1D("h4final","title;xaxis title;yaxis title",20,0,3.4);



    Long64_t nbytes = 0, nb = 0;
    
    
    std::cout << nentries << std::endl;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
       
       bool signal = (phoSigmaIEtaIEta_2012[0]<0.010);
       bool sideband = (phoSigmaIEtaIEta_2012[0]>0.011 && phoSigmaIEtaIEta_2012[0]<0.017);
       
       if (jentry % 1000 == 0) {
           std::cout << jentry << std::endl;}
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
       
       
       // select centrality 0-30% (hiBin between 0 and 60)
       // updated: select phoEtCorrected
       if ((hiBin < 0 )|| (hiBin > 60)) {
           continue;}
       /* // probably need to loop through jet arrays and individually select jets
       for (int jetloop = 0; jetloop < 9; jetloop++):
           if ;*/
       
       for (Long64_t jetTrk = 0; jetTrk < njet; jetTrk++) {
           if ((jetpt[jetTrk] <= 30) || (std::abs(jeteta[jetTrk]) >= 1.6) || (std::abs(jetphi[jetTrk] - phoPhi[0]) <= (7*M_PI/8))){
               continue;}
           
           //put into different plots based on centrality
           
           if ((hiBin >= 50.0) && (hiBin < 100.0)) {
               h1->Fill(jetpt[jetTrk] / phoEtCorrected[0]);
               h1m->Fill(jetpt_mix[jetTrk] / phoEtCorrected[0]);}
           if ((hiBin >= 30.0) && (hiBin < 50.0)) {
               h2->Fill(jetpt[jetTrk] / phoEtCorrected[0]);
               h2m->Fill(jetpt_mix[jetTrk] / phoEtCorrected[0]);}
           if ((hiBin >= 10.0) && (hiBin < 30.0)) {
               h3->Fill(jetpt[jetTrk] / phoEtCorrected[0]);
               h3m->Fill(jetpt_mix[jetTrk] / phoEtCorrected[0]);}
           if ((hiBin > 0) && (hiBin < 10.00)) {
               h4->Fill(jetpt[jetTrk] / phoEtCorrected[0]);
               h4m->Fill(jetpt_mix[jetTrk] / phoEtCorrected[0]);
               /*std::cout << "h4 test" << std::endl;
               std::cout << jetpt[jetTrk] / phoEtCorrected[0] << std::endl */ }
           
       }
   }
    
    // normalize the normal histograms
    Double_t scale1 = 1/h1->GetEntries();
    Double_t scale2 = 1/h2->GetEntries();
    Double_t scale3 = 1/h3->GetEntries();
    Double_t scale4 = 1/h4->GetEntries();

    std::cout << scale1 << std::endl;
    std::cout << h4->GetEntries() << std::endl;
    h1->Scale(scale1);
    h2->Scale(scale2);
    h3->Scale(scale3);
    h4->Scale(scale4);
    
    // normalize the mixed histograms, then scale by nmix
    
    Double_t scale1m = 1/h1m->GetEntries();
    Double_t scale2m = 1/h2m->GetEntries();
    Double_t scale3m = 1/h3m->GetEntries();
    Double_t scale4m = 1/h4m->GetEntries();
    
    h1m ->Scale(scale1m/nmix);
    h2m ->Scale(scale2m/nmix);
    h3m ->Scale(scale3m/nmix);
    h4m ->Scale(scale4m/nmix);
    
    // subtract the histograms from each other
    h1->Add(h1m,-1);
    h2->Add(h2m,-1);
    h3->Add(h3m,-1);
    h4->Add(h4m,-1);
    
    //find signal and background distributions
    
    //combine into final distributions



    
    
    outputfile->Write();
    outputfile->Close();
}

int main(int argc, char *argv[])
{
    analysis * ana = new analysis();
    ana->Loop();
    return 0;
}