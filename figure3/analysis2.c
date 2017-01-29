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
//    Long64_t ncentbins = 16;

    TFile * outputfile = new TFile("meantrkPt.root","recreate");
    /*
    TH1D * h1 = new TH1D("h1test","title;xaxis title;yaxis title",20,0,3.4);
    TH1D * h1m = new TH1D("h1mixed","title;xaxis title;yaxis title",20,0,3.4);
    TH1D * h2 = new TH1D("h2test","title;xaxis title;yaxis title",20,0,3);
    TH1D * h2m = new TH1D("h2mixed","title;xaxis title;yaxis title",20,0,3);
    TH1D * h3 = new TH1D("h3test","title;xaxis title;yaxis title",16,0,2);
    TH1D * h3m = new TH1D("h3mixed","title;xaxis title;yaxis title",16,0,2);
    TH1D * h4 = new TH1D("h4test","title;xaxis title;yaxis title",20,0,3.4);
    TH1D * h4m = new TH1D("h4mixed","title;xaxis title;yaxis title",20,0,3.4); */
    
    TH1D * h1sig = new TH1D("h1sig","title;xaxis title;yaxis title",16,0,2);
    TH1D * h1bg = new TH1D("h1bg","title;xaxis title;yaxis title",16,0,2);
    TH1D * h2sig = new TH1D("h2sig","title;xaxis title;yaxis title",16,0,2);
    TH1D * h2bg = new TH1D("h2bg","title;xaxis title;yaxis title",16,0,2);
    TH1D * h3sig = new TH1D("h3sig","title;xaxis title;yaxis title",16,0,2);
    TH1D * h3bg = new TH1D("h3bg","title;xaxis title;yaxis title",16,0,2);
    TH1D * h4sig = new TH1D("h4sig","title;xaxis title;yaxis title",16,0,2);
    TH1D * h4bg = new TH1D("h4bg","title;xaxis title;yaxis title",16,0,2);
    
    TH1D * h1msig = new TH1D("h1sig (mixed)","title;xaxis title;yaxis title",16,0,2);
    TH1D * h1mbg = new TH1D("h1bg (mixed)","title;xaxis title;yaxis title",16,0,2);
    TH1D * h2msig = new TH1D("h2sig (mixed)","title;xaxis title;yaxis title",16,0,2);
    TH1D * h2mbg = new TH1D("h2bg (mixed)","title;xaxis title;yaxis title",16,0,2);
    TH1D * h3msig = new TH1D("h3sig (mixed)","title;xaxis title;yaxis title",16,0,2);
    TH1D * h3mbg = new TH1D("h3bg (mixed)","title;xaxis title;yaxis title",16,0,2);
    TH1D * h4msig = new TH1D("h4sig (mixed)","title;xaxis title;yaxis title",16,0,2);
    TH1D * h4mbg = new TH1D("h4bg (mixed)","title;xaxis title;yaxis title",16,0,2);
    
    TH1D * h1f = new TH1D("h1final","title;xaxis title;yaxis title",16,0,2);
    TH1D * h2f = new TH1D("h2final","title;xaxis title;yaxis title",16,0,2);
    TH1D * h3f = new TH1D("h3final","title;xaxis title;yaxis title",16,0,2);
    TH1D * h4f = new TH1D("h4final","title;xaxis title;yaxis title",16,0,2);
    
    TH1D * phoEtSignal = new TH1D("phoEtSignal","title;xaxis title;yaxis title",16,0,200);
    TH1D * phoEtBackground = new TH1D("phoEtBg","title;xaxis title;yaxis title",16,0,200);


    Long64_t nmixTot = 0;
    float binwidth = h1f->GetBinWidth(1);
    
    Long64_t nbytes = 0, nb = 0;
    std::cout << nentries << std::endl;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        nmixTot += nmix;
        //std::cout << "nmix" << nmix << std::endl;
        if (jentry % 1000 == 0) {
            std::cout << jentry << std::endl;}
        if(jentry>30000) break;
        Long64_t ientry = LoadTree(jentry);
//        std::cout << "tree loaded" << std::endl;
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
       // std::cout<<nmix<<std::endl;
        bool signal = (phoSigmaIEtaIEta_2012[0]<0.010);
        bool sideband = (phoSigmaIEtaIEta_2012[0]>0.011 && phoSigmaIEtaIEta_2012[0]<0.017);
        // select signal
        
        if (signal) {
            //            std::cout << "sideband detected" << std::endl;
            // select centrality 0-30% (hiBin between 0 and 60)
            // updated: select phoEtCorrected
            //std::cout << "phoEtCorrected[0]" << phoEtCorrected[0] <<std::endl;
            //phoEtSignal->Fill(1);
            phoEtSignal->Fill(phoEtCorrected[0]*1.0);
            
            if (phoEtCorrected[0] < 60) {
                //                std::cout << "cutting phoEtCorrected 2" << std::endl;
                continue;}
            
            //            std::cout << "cut phoEtCorrected 2" << std::endl;
            // loop through jets
            //            std::cout << hiBin << std::endl;
            for (Long64_t jetTrk = 0; jetTrk < njet; jetTrk++) {
                if ((jetpt[jetTrk] <= 30) || (std::abs(jeteta[jetTrk]) >= 1.6) || (std::abs(jetphi[jetTrk] - phoPhi[0]) <= (7*M_PI/8))){
                    continue;}
                
                phoEtSignal->Fill(phoEtCorrected[0]);
                //                std::cout << "passed next cut 2" << std::endl;
                
                //put into different plots based on centrality
                
                if ((hiBin >= 100.0) && (hiBin < 200.0)) {
                    h1sig->Fill(jetpt[jetTrk] / phoEtCorrected[0]);
                }
                if ((hiBin >= 60.0) && (hiBin < 100.0)) {
                    h2sig->Fill(jetpt[jetTrk] / phoEtCorrected[0]);
                }
                if ((hiBin >= 20.0) && (hiBin < 60.0)) {
                    h3sig->Fill(jetpt[jetTrk] / phoEtCorrected[0]);
                }
                if ((hiBin > 0) && (hiBin < 20.00)) {
                    h4sig->Fill(jetpt[jetTrk] / phoEtCorrected[0]);
                }
                /*std::cout << "h4 test" << std::endl;
                 std::cout << jetpt[jetTrk] / phoEtCorrected[0] << std::endl */}
            
            // mixed jets
            
            for (Long64_t jetTrkMix = 0; jetTrkMix < njet_mix; jetTrkMix++) {
                if ((jetpt_mix[jetTrkMix] <= 30) || (std::abs(jeteta_mix[jetTrkMix]) >= 1.6) || (std::abs(jetphi_mix[jetTrkMix] - phoPhi[0]) <= (7*M_PI/8))){
                    continue;}
                if (nmix == 0) {
                    continue;
                }
                if ((hiBin >= 100.0) && (hiBin < 200.0)) {
                    h1msig->Fill(jetpt_mix[jetTrkMix] / (phoEtCorrected[0]),1/nmix);}
                if ((hiBin >= 60.0) && (hiBin < 100.0)) {
                    h2msig->Fill(jetpt_mix[jetTrkMix] / (phoEtCorrected[0]),1/nmix);}
                if ((hiBin >= 20.0) && (hiBin < 60.0)) {
                    h3msig->Fill(jetpt_mix[jetTrkMix] / (phoEtCorrected[0]),1/nmix);}
                if ((hiBin > 0) && (hiBin < 20.00)) {
                    h4msig->Fill(jetpt_mix[jetTrkMix] / (phoEtCorrected[0]),1/nmix);}
                
            }
        }
        
        if (sideband) {
//            std::cout << "sideband detected" << std::endl;
            // select centrality 0-30% (hiBin between 0 and 60)
            // updated: select phoEtCorrected
            //phoEtBackground->Fill(1);
            phoEtBackground->Fill(phoEtCorrected[0]*1.0);
            if (phoEtCorrected[0] < 60) {
//                std::cout << "cutting phoEtCorrected 2" << std::endl;
                continue;}
//            std::cout << "cut phoEtCorrected 2" << std::endl;
            // loop through jets
//            std::cout << hiBin << std::endl;
            for (Long64_t jetTrk = 0; jetTrk < njet; jetTrk++) {
                if ((jetpt[jetTrk] <= 30) || (std::abs(jeteta[jetTrk]) >= 1.6) || (std::abs(jetphi[jetTrk] - phoPhi[0]) <= (7*M_PI/8))){
                    continue;}
                
                phoEtBackground->Fill(phoEtCorrected[0]);

//                std::cout << "passed next cut 2" << std::endl;
                
                //put into different plots based on centrality
                
                if ((hiBin >= 100.0) && (hiBin < 200.0)) {
                    h1bg->Fill(jetpt[jetTrk] / phoEtCorrected[0]);
                }
                if ((hiBin >= 60.0) && (hiBin < 100.0)) {
                    h2bg->Fill(jetpt[jetTrk] / phoEtCorrected[0]);
                }
                if ((hiBin >= 20.0) && (hiBin < 60.0)) {
                    h3bg->Fill(jetpt[jetTrk] / phoEtCorrected[0]);
                }
                if ((hiBin > 0) && (hiBin < 20.00)) {
                    h4bg->Fill(jetpt[jetTrk] / phoEtCorrected[0]);
                }
                    /*std::cout << "h4 test" << std::endl;
                     std::cout << jetpt[jetTrk] / phoEtCorrected[0] << std::endl */}
            
            // mixed jets
            for (Long64_t jetTrkMix = 0; jetTrkMix < njet_mix; jetTrkMix++) {
//                if ((jetpt_mix[jetTrkMix] <= 30) || (std::abs(jeteta_mix[jetTrkMix]) >= 1.6) || (std::abs(jetphi_mix[jetTrkMix] - phoPhi[0]) <= (7*M_PI/8))){
//                    continue;}
                if (nmix == 0) {
                    continue;
                }
                if ((hiBin >= 100.0) && (hiBin < 200.0)) {
                    h1mbg->Fill(jetpt_mix[jetTrkMix] / (phoEtCorrected[0]),1/nmix);}
                if ((hiBin >= 60.0) && (hiBin < 100.0)) {
                    h2mbg->Fill(jetpt_mix[jetTrkMix] / (phoEtCorrected[0]),1/nmix);}
                if ((hiBin >= 20.0) && (hiBin < 60.0)) {
                    h3mbg->Fill(jetpt_mix[jetTrkMix] / (phoEtCorrected[0]),1/nmix);}
                if ((hiBin > 0) && (hiBin < 20.00)) {
                    h4mbg->Fill(jetpt_mix[jetTrkMix] / (phoEtCorrected[0]),1/nmix);}

            }
        }
    }
    
        // weight entries of mixed distribution by nmix
    // normalize the normal histograms
    
    std::cout << "h1sig before normalization" << h1sig->GetBinContent(2) << std::endl;
    std::cout << "binwidth" << binwidth << std::endl;
    binwidth = 0.125;
    

    
    std::cout << "h1sig after normalization" << h1sig->GetBinContent(2) << std::endl;

    
    // normalize the mixed histograms, then scale by nmix
    
    std::cout << "h1msig after scaling" << h1msig->GetBinContent(2) << std::endl;

    // subtract the histograms from each other
    h1sig->Add(h1msig,-1);
    h2sig->Add(h2msig,-1);
    h3sig->Add(h3msig,-1);
    h4sig->Add(h4msig,-1);
        
    h1bg->Add(h1mbg,-1);
    h2bg->Add(h2mbg,-1);
    h3bg->Add(h3mbg,-1);
    h4bg->Add(h4mbg,-1);
    
    
    std::cout << "h1sig after combination with h1msig" << h1sig->GetBinContent(2) << std::endl;

    //calculate nphosig, nphobg, binwidth, scale each histogram by 1/(npho*binwidth)
    float nphosig = phoEtSignal->Integral();
    float nphobg = phoEtBackground->Integral();
    
    std::cout << "binwidth" << binwidth << std::endl;
    std::cout << "nphosig" << nphosig << std::endl;
    std::cout << "nphobg" << nphobg << std::endl;
    
    
    h1sig->Scale(1/(binwidth*nphosig));
    h2sig->Scale(1/(binwidth*nphosig));
    h3sig->Scale(1/(binwidth*nphosig));
    h4sig->Scale(1/(binwidth*nphosig));
    
    h1bg->Scale(1/(binwidth * nphobg));
    h2bg->Scale(1/(binwidth * nphobg));
    h3bg->Scale(1/(binwidth * nphobg));
    h4bg->Scale(1/(binwidth * nphobg));
    
    //combine into final distributions (final = (1/purity) * signal - ((1-purity)/purity) * background)
    
    // calculate purity CHECK LAST ENTRY (ispp value)
        
        Double_t purity1 = getpurity(60,100,isPP);
        Double_t purity2 = getpurity(60,60,isPP);
        Double_t purity3 = getpurity(60,20,isPP);
        Double_t purity4 = getpurity(60,0,isPP);
    
        // test:
    std::cout << "purity1" << purity1 << std::endl;
    std::cout << "purity2" << purity2 << std::endl;
    std::cout << "purity3" << purity3 << std::endl;
    std::cout << "purity4" << purity4 << std::endl;
    
        h1sig->Scale(1/purity1);
        h1bg ->Scale(((1-purity1)/purity1));
        h1f->Add(h1sig,1);
        h1f->Add(h1bg,-1);
    
    std::cout << "h1sig after purity scaling" << h1sig->GetBinContent(2) << std::endl;

        h2sig->Scale(1/purity2);
        h2bg ->Scale(((1-purity2)/purity2));
        h2f->Add(h2sig,1);
        h2f->Add(h2bg,-1);

        h3sig->Scale(1/purity3);
        h3bg ->Scale(((1-purity3)/purity3));
        h3f->Add(h3sig,1);
        h3f->Add(h3bg,-1);

        h4sig->Scale(1/purity4);
        h4bg ->Scale(((1-purity4)/purity4));
        h4f->Add(h4sig,1);
        h4f->Add(h4bg,-1);
    
////     normalize final histograms
//    Double_t scale1f = 1/h1f->GetEntries();
//    Double_t scale2f = 1/h2f->GetEntries();
//    Double_t scale3f = 1/h3f->GetEntries();
//    Double_t scale4f = 1/h4f->GetEntries();
//    
////     instead of scaling by total number of events,
//    h1f ->Scale(scale1f);
//    h2f ->Scale(scale2f);
//    h3f ->Scale(scale3f);
//    h4f ->Scale(scale4f);
    
    std::cout << h1f->GetBinContent(2) << std::endl;
    
    
    //std::string cents[] = {"0-10%","10-30%","30-50%","50-100%"};
    
    outputfile->Write();
    outputfile->Close();

}

int main(int argc, char *argv[])
{
    analysis * ana = new analysis();
    ana->Loop();
    return 0;
}