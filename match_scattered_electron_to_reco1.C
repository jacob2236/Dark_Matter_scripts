/*
 * match_scattered_electron_to_reco.C
 *
 * Matches the MC scattered electron to its reconstructed charged particle
 * in reco.root (EICrecon output).
 *
 * MATCHING CHAIN:
 *  1. MCScatteredElectrons_objIdx.index
 *       -> MCParticle index of the EICrecon-identified scattered electron
 *  2. Search _ReconstructedChargedParticleAssociations_sim.index for that MC index
 *       -> gives the row in the association table
 *  3. Read _ReconstructedChargedParticleAssociations_rec.index at that row
 *       -> index into ReconstructedChargedParticles
 *  4. Look up ReconstructedChargedParticles at that index
 *       -> reco momentum, energy, charge, PDG, mass
 *
 *  If multiple association entries match (shared MC particle), pick the one
 *  with the highest weight.
 *
 *  Fallback: if MCScatteredElectrons_objIdx is empty, fall back to the
 *  original status-23 + status-21-parent method from MCParticles.
 *
 * Output:
 *   - reco_scattered_electron_output.root
 *   - reco_scattered_electron_kinematics.png/.pdf
 *
 * Usage: root -l -b -q 'match_scattered_electron_to_reco.C("reco.root")'
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"

using namespace std;

// -----------------------------------------------------------------------
// Safe leaf getter
// -----------------------------------------------------------------------
TLeaf* getLeaf(TTree* tree, const char* name, bool required = true) {
    TLeaf* l = tree->GetLeaf(name);
    if (!l && required) cerr << "  ERROR: required leaf missing: " << name << endl;
    else if (!l)        cerr << "  WARNING: optional leaf missing: " << name << endl;
    return l;
}

// -----------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------
void match_scattered_electron_to_reco1(const char* filename = "reco.root") {

    cout << "\n============================================================" << endl;
    cout << "  MC SCATTERED ELECTRON -> RECO PARTICLE MATCHER" << endl;
    cout << "============================================================\n" << endl;

    TFile* f = TFile::Open(filename);
    if (!f || f->IsZombie()) { cerr << "Cannot open: " << filename << endl; return; }
    TTree* tree = (TTree*)f->Get("events");
    if (!tree) { cerr << "No 'events' tree." << endl; return; }

    Long64_t nEvents = tree->GetEntries();
    cout << "File  : " << filename << endl;
    cout << "Events: " << nEvents  << "\n" << endl;

    // ----------------------------------------------------------------
    // MCParticles leaves  (for fallback + MC truth kinematics)
    // ----------------------------------------------------------------
    TLeaf* mc_pdg    = getLeaf(tree, "MCParticles.PDG");
    TLeaf* mc_status = getLeaf(tree, "MCParticles.generatorStatus");
    TLeaf* mc_px     = getLeaf(tree, "MCParticles.momentum.x");
    TLeaf* mc_py     = getLeaf(tree, "MCParticles.momentum.y");
    TLeaf* mc_pz     = getLeaf(tree, "MCParticles.momentum.z");
    TLeaf* mc_mass   = getLeaf(tree, "MCParticles.mass");
    TLeaf* mc_pb     = getLeaf(tree, "MCParticles.parents_begin");
    TLeaf* mc_pe     = getLeaf(tree, "MCParticles.parents_end");
    TLeaf* mc_paridx = getLeaf(tree, "_MCParticles_parents.index", false);

    // ----------------------------------------------------------------
    // EICrecon scattered electron identification
    // MCScatteredElectrons_objIdx.index  -> MCParticle index
    // ----------------------------------------------------------------
    TLeaf* scEl_mcIdx = getLeaf(tree, "MCScatteredElectrons_objIdx.index", false);

    // ----------------------------------------------------------------
    // MC <-> Reco association
    // _ReconstructedChargedParticleAssociations_sim.index -> MC index
    // _ReconstructedChargedParticleAssociations_rec.index -> Reco index
    // ReconstructedChargedParticleAssociations.weight     -> match weight
    // ----------------------------------------------------------------
    TLeaf* assoc_sim    = getLeaf(tree, "_ReconstructedChargedParticleAssociations_sim.index");
    TLeaf* assoc_rec    = getLeaf(tree, "_ReconstructedChargedParticleAssociations_rec.index");
    TLeaf* assoc_weight = getLeaf(tree, "ReconstructedChargedParticleAssociations.weight", false);

    // ----------------------------------------------------------------
    // ReconstructedChargedParticles leaves
    // ----------------------------------------------------------------
    TLeaf* reco_px     = getLeaf(tree, "ReconstructedChargedParticles.momentum.x");
    TLeaf* reco_py     = getLeaf(tree, "ReconstructedChargedParticles.momentum.y");
    TLeaf* reco_pz     = getLeaf(tree, "ReconstructedChargedParticles.momentum.z");
    TLeaf* reco_energy = getLeaf(tree, "ReconstructedChargedParticles.energy");
    TLeaf* reco_charge = getLeaf(tree, "ReconstructedChargedParticles.charge");
    TLeaf* reco_mass   = getLeaf(tree, "ReconstructedChargedParticles.mass");
    TLeaf* reco_pdg    = getLeaf(tree, "ReconstructedChargedParticles.PDG", false);

    if (!mc_pdg || !assoc_sim || !assoc_rec || !reco_px || !reco_py || !reco_pz || !reco_energy) {
        cerr << "Missing required leaves — cannot continue." << endl;
        return;
    }

    // ----------------------------------------------------------------
    // Output
    // ----------------------------------------------------------------
    TFile* outFile = new TFile("reco_scattered_electron_output.root", "RECREATE");
    TTree* outTree = new TTree("scattered_electrons", "MC + Reco scattered electron");

    int    evtNum;
    // MC truth
    int    mc_idx, mc_pdg_b;
    double mc_pt, mc_eta, mc_phi, mc_energy, mc_p;
    bool   mc_found, mc_used_eicrecon, mc_used_fallback;
    // Reco match
    int    reco_idx, reco_pdg_b;
    double reco_pt_b, reco_eta_b, reco_phi_b, reco_energy_b, reco_p_b;
    float  reco_charge_b, reco_mass_b;
    bool   reco_found;
    double assoc_w;
    // Residuals
    double delta_pt, delta_eta, delta_phi, delta_p, delta_E;

    outTree->Branch("event",           &evtNum,           "event/I");
    // MC
    outTree->Branch("mc_found",        &mc_found,         "mc_found/O");
    outTree->Branch("mc_eicrecon",     &mc_used_eicrecon, "mc_eicrecon/O");
    outTree->Branch("mc_fallback",     &mc_used_fallback, "mc_fallback/O");
    outTree->Branch("mc_idx",          &mc_idx,           "mc_idx/I");
    outTree->Branch("mc_pdg",          &mc_pdg_b,         "mc_pdg/I");
    outTree->Branch("mc_pt",           &mc_pt,            "mc_pt/D");
    outTree->Branch("mc_eta",          &mc_eta,           "mc_eta/D");
    outTree->Branch("mc_phi",          &mc_phi,           "mc_phi/D");
    outTree->Branch("mc_energy",       &mc_energy,        "mc_energy/D");
    outTree->Branch("mc_p",            &mc_p,             "mc_p/D");
    // Reco
    outTree->Branch("reco_found",      &reco_found,       "reco_found/O");
    outTree->Branch("reco_idx",        &reco_idx,         "reco_idx/I");
    outTree->Branch("reco_pdg",        &reco_pdg_b,       "reco_pdg/I");
    outTree->Branch("reco_pt",         &reco_pt_b,        "reco_pt/D");
    outTree->Branch("reco_eta",        &reco_eta_b,       "reco_eta/D");
    outTree->Branch("reco_phi",        &reco_phi_b,       "reco_phi/D");
    outTree->Branch("reco_energy",     &reco_energy_b,    "reco_energy/D");
    outTree->Branch("reco_p",          &reco_p_b,         "reco_p/D");
    outTree->Branch("reco_charge",     &reco_charge_b,    "reco_charge/F");
    outTree->Branch("reco_mass",       &reco_mass_b,      "reco_mass/F");
    outTree->Branch("assoc_weight",    &assoc_w,          "assoc_weight/D");
    // Residuals (reco - MC)
    outTree->Branch("delta_pt",        &delta_pt,         "delta_pt/D");
    outTree->Branch("delta_eta",       &delta_eta,        "delta_eta/D");
    outTree->Branch("delta_phi",       &delta_phi,        "delta_phi/D");
    outTree->Branch("delta_p",         &delta_p,          "delta_p/D");
    outTree->Branch("delta_E",         &delta_E,          "delta_E/D");

    // ----------------------------------------------------------------
    // Histograms
    // ----------------------------------------------------------------
    // MC truth
    TH1F* h_mc_pt     = new TH1F("h_mc_pt",    "MC Scattered e^{-} p_{T};p_{T} (GeV);Events",    60, 0, 25);
    TH1F* h_mc_eta    = new TH1F("h_mc_eta",   "MC Scattered e^{-} #eta;#eta;Events",             60,-5, 5);
    TH1F* h_mc_energy = new TH1F("h_mc_energy","MC Scattered e^{-} Energy;E (GeV);Events",        60, 0, 30);
    // Reco
    TH1F* h_reco_pt     = new TH1F("h_reco_pt",    "Reco Scattered e^{-} p_{T};p_{T} (GeV);Events",   60, 0, 25);
    TH1F* h_reco_eta    = new TH1F("h_reco_eta",   "Reco Scattered e^{-} #eta;#eta;Events",            60,-5, 5);
    TH1F* h_reco_energy = new TH1F("h_reco_energy","Reco Scattered e^{-} Energy;E (GeV);Events",       60, 0, 30);
    // Residuals
    TH1F* h_dpt   = new TH1F("h_dpt",  "Reco - MC p_{T};#Deltap_{T} (GeV);Events",    80,-4, 4);
    TH1F* h_deta  = new TH1F("h_deta", "Reco - MC #eta;#Delta#eta;Events",             80,-0.5,0.5);
    TH1F* h_dphi  = new TH1F("h_dphi", "Reco - MC #phi;#Delta#phi (rad);Events",       80,-0.5,0.5);
    TH1F* h_dE    = new TH1F("h_dE",   "Reco - MC Energy;#DeltaE (GeV);Events",        80,-5, 5);
    TH1F* h_dp    = new TH1F("h_dp",   "Reco - MC |p|;#Delta|p| (GeV);Events",         80,-5, 5);
    // Relative residuals
    TH1F* h_dpt_rel  = new TH1F("h_dpt_rel",  "(Reco-MC)/MC p_{T};(p_{T}^{reco}-p_{T}^{MC})/p_{T}^{MC};Events", 80,-0.5,0.5);
    TH1F* h_dE_rel   = new TH1F("h_dE_rel",   "(Reco-MC)/MC Energy;(E^{reco}-E^{MC})/E^{MC};Events",             80,-0.5,0.5);
    // 2D: reco vs MC
    TH2F* h_pt2d  = new TH2F("h_pt2d",  "Reco p_{T} vs MC p_{T};p_{T}^{MC} (GeV);p_{T}^{reco} (GeV)",
                               60, 0,25, 60, 0,25);
    TH2F* h_eta2d = new TH2F("h_eta2d", "Reco #eta vs MC #eta;#eta^{MC};#eta^{reco}",
                               60,-5,5, 60,-5,5);
    TH2F* h_E2d   = new TH2F("h_E2d",   "Reco E vs MC E;E^{MC} (GeV);E^{reco} (GeV)",
                               60, 0,30, 60, 0,30);
    // Matching weight & reco charge
    TH1F* h_weight = new TH1F("h_weight", "Association weight;weight;Entries", 50, 0, 1.1);
    TH1F* h_charge = new TH1F("h_charge", "Reco charge;charge;Entries",        7,-3.5,3.5);

    // ----------------------------------------------------------------
    // Counters
    // ----------------------------------------------------------------
    int nMCFound=0, nRecoFound=0, nEICRecon=0, nFallback=0;

    // ================================================================
    // EVENT LOOP
    // ================================================================
    cout << "Processing events..." << endl;
    cout << "------------------------------------------------------------" << endl;

    for (Long64_t iEvt = 0; iEvt < nEvents; ++iEvt) {
        tree->GetEntry(iEvt);

        evtNum          = (int)iEvt;
        mc_found        = false;
        mc_used_eicrecon= false;
        mc_used_fallback= false;
        reco_found      = false;
        mc_idx          = -1;  mc_pdg_b = 0;
        mc_pt = mc_eta = mc_phi = mc_energy = mc_p = 0;
        reco_idx        = -1;  reco_pdg_b = 0;
        reco_pt_b = reco_eta_b = reco_phi_b = reco_energy_b = reco_p_b = 0;
        reco_charge_b   = 0;   reco_mass_b = 0;
        assoc_w         = 0;
        delta_pt = delta_eta = delta_phi = delta_p = delta_E = 0;

        int nParticles = mc_pdg->GetLen();

        // --------------------------------------------------------
        // STEP 1: Get MC scattered electron index
        //
        // KEY FINDING from debug: MCScatteredElectrons_objIdx always has
        // TWO entries — one PDG=11 (scattered electron) and one PDG=-11
        // (beam positron remnant). We must select PDG=11 only.
        //
        // The old status-23+parent fallback found a DIFFERENT MC index
        // (the pre-scatter beam electron) which never appears in the
        // association table. Always use MCScatteredElectrons_objIdx.
        // --------------------------------------------------------
        int selMCIdx = -1;

        // Method A: EICrecon MCScatteredElectrons — pick the PDG=11 entry
        if (scEl_mcIdx && scEl_mcIdx->GetLen() > 0) {
            double maxE = -1;
            for (int i = 0; i < scEl_mcIdx->GetLen(); ++i) {
                int idx = (int)scEl_mcIdx->GetValue(i);
                if (idx < 0 || idx >= nParticles) continue;
                // Must be electron (PDG=11), not beam positron (PDG=-11)
                if ((int)mc_pdg->GetValue(idx) != 11) continue;
                double px = mc_px->GetValue(idx);
                double py = mc_py->GetValue(idx);
                double pz = mc_pz->GetValue(idx);
                double m  = mc_mass->GetValue(idx);
                double e  = sqrt(px*px+py*py+pz*pz+m*m);
                if (e > maxE) { maxE = e; selMCIdx = idx; }
            }
            if (selMCIdx >= 0) {
                mc_used_eicrecon = true;
                nEICRecon++;
            }
        }

        // Method B: fallback — if MCScatteredElectrons_objIdx gave nothing,
        // find any PDG=11 status-23 particle that appears in the association
        // table (guarantees the index is valid for reco matching).
        if (selMCIdx < 0) {
            int nAssocFB = assoc_sim->GetLen();
            double maxE = -1;
            for (int i = 0; i < nParticles; ++i) {
                if ((int)mc_pdg->GetValue(i) != 11)   continue;
                if ((int)mc_status->GetValue(i) != 23) continue;
                // Only consider if this index appears in the association table
                for (int a = 0; a < nAssocFB; ++a) {
                    if ((int)assoc_sim->GetValue(a) == i) {
                        double px = mc_px->GetValue(i);
                        double py = mc_py->GetValue(i);
                        double pz = mc_pz->GetValue(i);
                        double m  = mc_mass->GetValue(i);
                        double e  = sqrt(px*px+py*py+pz*pz+m*m);
                        if (e > maxE) { maxE = e; selMCIdx = i; }
                        break;
                    }
                }
            }
            if (selMCIdx >= 0) {
                mc_used_fallback = true;
                nFallback++;
            }
        }

        if (selMCIdx < 0) { outTree->Fill(); continue; }

        // --------------------------------------------------------
        // Fill MC truth kinematics
        // --------------------------------------------------------
        {
            double px = mc_px->GetValue(selMCIdx);
            double py = mc_py->GetValue(selMCIdx);
            double pz = mc_pz->GetValue(selMCIdx);
            double m  = mc_mass->GetValue(selMCIdx);
            double p  = sqrt(px*px+py*py+pz*pz);
            double e  = sqrt(p*p+m*m);
            double pt = sqrt(px*px+py*py);
            double eta = (p>0 && p!=fabs(pz)) ? 0.5*log((p+pz)/(p-pz)) : 0;
            double phi = atan2(py,px);

            mc_found   = true;
            mc_idx     = selMCIdx;
            mc_pdg_b   = (int)mc_pdg->GetValue(selMCIdx);
            mc_pt      = pt;
            mc_eta     = eta;
            mc_phi     = phi;
            mc_energy  = e;
            mc_p       = p;

            h_mc_pt->Fill(pt);
            h_mc_eta->Fill(eta);
            h_mc_energy->Fill(e);
            nMCFound++;
        }

        // --------------------------------------------------------
        // STEP 2: Find matching reco particle via association table
        // --------------------------------------------------------
        int    bestRecoIdx = -1;
        double bestWeight  = -1;
        int    nAssoc      = assoc_sim->GetLen();

        for (int a = 0; a < nAssoc; ++a) {
            if ((int)assoc_sim->GetValue(a) != selMCIdx) continue;
            double w = assoc_weight ? (double)assoc_weight->GetValue(a) : 1.0;
            if (w > bestWeight) {
                bestWeight  = w;
                bestRecoIdx = (int)assoc_rec->GetValue(a);
            }
        }

        if (bestRecoIdx < 0 || bestRecoIdx >= reco_px->GetLen()) {
            outTree->Fill();
            continue;
        }

        // --------------------------------------------------------
        // Fill reco kinematics
        // --------------------------------------------------------
        {
            double px  = reco_px->GetValue(bestRecoIdx);
            double py  = reco_py->GetValue(bestRecoIdx);
            double pz  = reco_pz->GetValue(bestRecoIdx);
            double e   = reco_energy->GetValue(bestRecoIdx);
            double p   = sqrt(px*px+py*py+pz*pz);
            double pt  = sqrt(px*px+py*py);
            double eta = (p>0 && p!=fabs(pz)) ? 0.5*log((p+pz)/(p-pz)) : 0;
            double phi = atan2(py,px);

            reco_found     = true;
            reco_idx       = bestRecoIdx;
            reco_pdg_b     = reco_pdg ? (int)reco_pdg->GetValue(bestRecoIdx) : 0;
            reco_pt_b      = pt;
            reco_eta_b     = eta;
            reco_phi_b     = phi;
            reco_energy_b  = e;
            reco_p_b       = p;
            reco_charge_b  = reco_charge ? (float)reco_charge->GetValue(bestRecoIdx) : 0;
            reco_mass_b    = reco_mass   ? (float)reco_mass->GetValue(bestRecoIdx)   : 0;
            assoc_w        = bestWeight;

            // Residuals
            delta_pt  = pt  - mc_pt;
            delta_eta = eta - mc_eta;
            delta_phi = phi - mc_phi;
            // Wrap delta_phi to [-pi, pi]
            while (delta_phi >  M_PI) delta_phi -= 2*M_PI;
            while (delta_phi < -M_PI) delta_phi += 2*M_PI;
            delta_p   = p   - mc_p;
            delta_E   = e   - mc_energy;

            h_reco_pt->Fill(pt);
            h_reco_eta->Fill(eta);
            h_reco_energy->Fill(e);
            h_dpt->Fill(delta_pt);
            h_deta->Fill(delta_eta);
            h_dphi->Fill(delta_phi);
            h_dE->Fill(delta_E);
            h_dp->Fill(delta_p);
            if (mc_pt    > 0) h_dpt_rel->Fill(delta_pt  / mc_pt);
            if (mc_energy> 0) h_dE_rel->Fill( delta_E   / mc_energy);
            h_pt2d->Fill(mc_pt,    pt);
            h_eta2d->Fill(mc_eta,  eta);
            h_E2d->Fill(mc_energy, e);
            h_weight->Fill(bestWeight);
            h_charge->Fill(reco_charge_b);

            nRecoFound++;
        }

        outTree->Fill();

        // Print first 10 matched events
        if (nRecoFound <= 10) {
            cout << "Event " << iEvt
                 << "  MC[" << selMCIdx << "]"
                 << (mc_used_eicrecon ? " [EICrecon]" : " [fallback]") << endl;
            cout << "  MC  : E=" << mc_energy << " GeV  pt=" << mc_pt
                 << "  eta=" << mc_eta << "  phi=" << mc_phi << endl;
            cout << "  Reco[" << bestRecoIdx << "]: E=" << reco_energy_b
                 << " GeV  pt=" << reco_pt_b
                 << "  eta=" << reco_eta_b << "  phi=" << reco_phi_b
                 << "  charge=" << reco_charge_b
                 << "  weight=" << bestWeight << endl;
            cout << "  dpt=" << delta_pt << "  deta=" << delta_eta
                 << "  dE=" << delta_E << endl;
        }

        if ((iEvt+1)%500==0 || iEvt==nEvents-1)
            cout << "  Processed " << (iEvt+1) << "/" << nEvents << " events\r" << flush;
    }

    // ================================================================
    // Summary
    // ================================================================
    cout << "\n\n============================================================" << endl;
    cout << "  RESULTS SUMMARY" << endl;
    cout << "============================================================" << endl;
    cout << "Total events          : " << nEvents   << endl;
    cout << "MC electron found     : " << nMCFound  << " (" << 100.0*nMCFound/nEvents  << "%)" << endl;
    cout << "  via EICrecon MCScatteredElectrons : " << nEICRecon  << endl;
    cout << "  via status-23 fallback            : " << nFallback  << endl;
    cout << "Reco match found      : " << nRecoFound<< " (" << (nMCFound>0?100.0*nRecoFound/nMCFound:0) << "% of MC found)" << endl;

    if (nRecoFound > 0) {
        cout << "\nMC truth kinematics:" << endl;
        cout << "  Mean pT    : " << h_mc_pt->GetMean()     << " +/- " << h_mc_pt->GetRMS()     << " GeV" << endl;
        cout << "  Mean eta   : " << h_mc_eta->GetMean()    << " +/- " << h_mc_eta->GetRMS()             << endl;
        cout << "  Mean Energy: " << h_mc_energy->GetMean() << " +/- " << h_mc_energy->GetRMS() << " GeV" << endl;
        cout << "\nReco kinematics:" << endl;
        cout << "  Mean pT    : " << h_reco_pt->GetMean()     << " +/- " << h_reco_pt->GetRMS()     << " GeV" << endl;
        cout << "  Mean eta   : " << h_reco_eta->GetMean()    << " +/- " << h_reco_eta->GetRMS()             << endl;
        cout << "  Mean Energy: " << h_reco_energy->GetMean() << " +/- " << h_reco_energy->GetRMS() << " GeV" << endl;
        cout << "\nResiduals (Reco - MC):" << endl;
        cout << "  Mean dpt   : " << h_dpt->GetMean()  << " +/- " << h_dpt->GetRMS()  << " GeV" << endl;
        cout << "  Mean deta  : " << h_deta->GetMean() << " +/- " << h_deta->GetRMS()          << endl;
        cout << "  Mean dE    : " << h_dE->GetMean()   << " +/- " << h_dE->GetRMS()   << " GeV" << endl;
        cout << "  Mean rel dpt: " << h_dpt_rel->GetMean() << " +/- " << h_dpt_rel->GetRMS()   << endl;
        cout << "  Mean rel dE : " << h_dE_rel->GetMean()  << " +/- " << h_dE_rel->GetRMS()    << endl;
        cout << "  Mean assoc weight: " << h_weight->GetMean() << endl;
    }
    cout << "============================================================\n" << endl;

    // ================================================================
    // Save
    // ================================================================
    outFile->cd();
    outTree->Write();

    // Canvas 1: kinematic comparison
    TCanvas* c1 = new TCanvas("c1","Kinematics Comparison",1800,1000);
    c1->Divide(3,2);
    auto styleStack = [](TH1F* hmc, TH1F* hreco) {
        hmc->SetLineColor(kBlue);   hmc->SetLineWidth(2);
        hreco->SetLineColor(kRed);  hreco->SetLineWidth(2);
        hreco->SetLineStyle(2);
    };
    c1->cd(1); styleStack(h_mc_pt,     h_reco_pt);     h_mc_pt->Draw();     h_reco_pt->Draw("SAME");
    c1->cd(2); styleStack(h_mc_eta,    h_reco_eta);    h_mc_eta->Draw();    h_reco_eta->Draw("SAME");
    c1->cd(3); styleStack(h_mc_energy, h_reco_energy); h_mc_energy->Draw(); h_reco_energy->Draw("SAME");
    c1->cd(4); h_pt2d->Draw("COLZ");
    c1->cd(5); h_eta2d->Draw("COLZ");
    c1->cd(6); h_E2d->Draw("COLZ");
    c1->Write();
    c1->SaveAs("reco_scattered_electron_kinematics.png");
    c1->SaveAs("reco_scattered_electron_kinematics.pdf");

    // Canvas 2: residuals
    TCanvas* c2 = new TCanvas("c2","Residuals",1800,600);
    c2->Divide(5,1);
    auto styleRes = [](TH1F* h){ h->SetLineColor(kBlack); h->SetLineWidth(2); h->SetFillColorAlpha(kAzure+1,0.3); };
    c2->cd(1); styleRes(h_dpt);  h_dpt->Draw();
    c2->cd(2); styleRes(h_deta); h_deta->Draw();
    c2->cd(3); styleRes(h_dphi); h_dphi->Draw();
    c2->cd(4); styleRes(h_dE);   h_dE->Draw();
    c2->cd(5); styleRes(h_dpt_rel); h_dpt_rel->Draw();
    c2->Write();
    c2->SaveAs("reco_scattered_electron_residuals.png");

    // Write all histograms
    for (auto h : {h_mc_pt, h_mc_eta, h_mc_energy,
                   h_reco_pt, h_reco_eta, h_reco_energy,
                   h_dpt, h_deta, h_dphi, h_dE, h_dp,
                   h_dpt_rel, h_dE_rel, h_weight, h_charge})
        h->Write();
    h_pt2d->Write(); h_eta2d->Write(); h_E2d->Write();

    outFile->Close();
    f->Close();

    cout << "Output files:" << endl;
    cout << "  reco_scattered_electron_output.root" << endl;
    cout << "  reco_scattered_electron_kinematics.png/.pdf" << endl;
    cout << "  reco_scattered_electron_residuals.png"       << endl;
    cout << "============================================================\n" << endl;
}
