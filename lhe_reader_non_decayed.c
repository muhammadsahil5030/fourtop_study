//&>/dev/null;x="${0%.*}";[ ! "$x" -ot "$0" ]||(rm -f "$x";g++ -o "$x" "$0" -I`root-config --incdir` `root-config --libs`);exit

// Build: g++ lhe_reader_non_decayed.c -o lhe_reader_non_decayed -I`root-config --incdir` `root-config --libs`
// OR, g++ -std=c++11 lhe_reader_non_decayed.c -o lhe_reader_non_decayed -I`root-config --incdir` `root-config --libs`
#include <cmath>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
//#include </home/muhammad/root6/build/include/TLorentzVector.h>
#include <TLorentzVector.h>
using namespace std;

// Pour une description du format leshouches
// hep-ph/0609017
//
// pour chaque evenement
// une ligne générale : NbPart idprocess poids scale alpha_em alpha_s
// pour chaque particule : id status mere1 mere2 couleur1 couleur2 px py pz E m lifetime spin  
int main(int argc, char **argv) {

  if (argc != 2) {
    cout << "Missing argument. Expected LHE filename without '.lhe'"<< endl;
    exit(-5);
  }
  float DeltaR(float eta1, float eta2, float phi1, float phi2);
  float DeltaPhi(float phi1, float phi2);

  string basename=argv[1];

  string lhefname = basename+".lhe";
  string rootfname = basename+".root";

  string tt;
  int event=0;
  int npart,idprocess;
  double weight,scale,alpha_em,alpha_s;

  TH1::SetDefaultSumw2(true);
//--------------------- insert here histogram definition -----------------------------------//

  TH1F* hPt_lep = new TH1F("lep_pt","Lepton pt",50,0.0,400) ;
  TH1F* hEta_lep = new TH1F("lep_eta","Lepton eta",50,-6,6) ;
  TH1F* hPhi_lep = new TH1F("lep_phi","Lepton phi",50,-3.16,3.16) ;
  TH1F* hPx_lep = new TH1F("lep_px","Lepton px",50,-200,200) ;
  TH1F* hPy_lep = new TH1F("lep_py","Lepton py",50,-200,200) ;

  TH1F* hPt_alep = new TH1F("alep_pt","anti-lepton pt",50,0.0,400) ;
  TH1F* hEta_alep = new TH1F("alep_eta","anti-lepton eta",50,-6,6) ;
  TH1F* hPhi_alep = new TH1F("alep_phi","anti-lepton phi",50,-3.16,3.16) ;
  TH1F* hPx_alep = new TH1F("alep_px","anti-lepton px",50,-200,200) ;
  TH1F* hPy_alep = new TH1F("alep_py","anti-lepton py",50,-200,200) ;
  
  TH1F* hPt_bJet = new TH1F("bJet_pt","Beauty Jet Pt",50,0.0,400) ;
  TH1F* hEta_bJet = new TH1F("bJet_eta","Beauty Jet Eta",50,-6,6) ;
  TH1F* hPhi_bJet = new TH1F("bJet_phi","Beauty Jet Phi",50,-3.16,3.16) ;
  
  TH1F* hPt_wboson = new TH1F("wboson_Pt","wboson Pt",50,0.0,500) ;
  TH1F* hEta_wboson = new TH1F("wboson_eta","wboson Eta",50,-6,6) ;
  TH1F* hPhi_wboson = new TH1F("wboson_Phi","wboson Phi",50,-3.16,3.16) ;
  
  TH1F* hM_wboson1 = new TH1F("wboson1_mass", "w1 boson mass", 50, 60, 100);
  TH1F* hM_wboson2 = new TH1F("wboson2_mass", "w2 boson mass", 50, 60, 100);
  TH1F* hM_wboson3 = new TH1F("wboson3_mass", "w3 boson mass", 50, 60, 100);
  TH1F* hM_wboson4 = new TH1F("wboson4_mass", "w4 boson mass", 50, 60, 100);
  TH1F* hM_wboson5 = new TH1F("wboson5_mass", "w5 boson mass", 50, 60, 100);
  
  TH1F* hM_top1 = new TH1F("top1_mass", "Top1 Quark mass", 50, 150, 200);
  TH1F* hM_top2 = new TH1F("top2_mass", "Top2 Quark mass", 50, 150, 200);
  TH1F* hM_top3 = new TH1F("top3_mass", "top3 Quark mass", 50, 150, 200);
  TH1F* hM_top4 = new TH1F("top4_mass", "top4 Quark mass", 50, 150, 200);
  TH1F* hM_top5 = new TH1F("top5_mass", "Top5 Quark mass", 50, 150, 200);

//------------------------------------------------------------------------------------------------------//  
// --- end histogram definition
cout<<"I am here: "<<endl;

  int nlept=0, nsemi=0, nhadr=0;
  ifstream ff(lhefname.c_str(),ios::in); //ouverture du fichier .lhe
  //ifstream ff("test.lhe",ios::in); //ouverture du fichier .lhe
  //ifstream ff("/home/cms/perries/madgraph-newphysics/s1/madevent/Events/zp4000_unweighted_events.lhe",ios::in); //ouverture du fichier .lhe
  //ifstream ff("/home/cms/perries/madgraph-newphysics/QCD/madevent/Events/qcd_unweighted_events.lhe",ios::in); //ouverture du fichier .lhe
  int negativeWeight = 0;
  long long line = 0;
  while(!ff.eof()) {
    std::stringstream buffer;
    ff>>tt;
    buffer << tt << std::endl;
    line++;

    if(tt=="<event>") {
      ff>>npart>>idprocess>>weight>>scale>>alpha_em>>alpha_s; //event definition
      buffer << npart << " " << idprocess << " " << weight << " " << scale << " " << alpha_em << " " << alpha_s << std::endl;
      line++;
      event++;
      if (weight < 0) {
        negativeWeight++;
        weight = -1;
      } else {
        weight = 1;
      }
      /*weight = 1.;*/

      if (event%1==0) cout << "reading event "<< event << endl;
      int lmin=-1, lmax=-1, bmin=-1, met1=-1, met2=-1, bmax=-1, met=-1, topLep,topbJet, tbarLep, tbarbJet, top_jet, tbar_jet;
      int Part_Id, Moth_Id, Part_absId, Moth_absId;
int n_lep=0, n_alep=0, n_lep_nu=0, n_topjets=0, n_tbarjets=0, n_topbjets=0, n_bJets=0, n_tbar_nu=0, n_top_nu=0, n_wlep=0;
int n_elec=0, n_elec_nu=0, n_muon=0, n_muon_nu=0;

int  bj[2]={-1,-1};
int lp[2]={-1,-1};

      int q[4]={-1,-1,-1,-1,};
      int muon[5]={-1,-1,-1,-1,-1};
      int elec[5]={-1,-1,-1,-1,-1};
      int elec_nu[5]={-1,-1,-1,-1,-1};
      int muon_nu[5]={-1,-1,-1,-1,-1};
      int bq[4]={-1,-1,-1,-1};
      
      int top=-1,topbar=-1,zprime=-1;
      int *Id      = new int[npart+1];
      int *Status  = new int[npart+1];
      int *Mother1 = new int[npart+1];
      int *Mother2 = new int[npart+1];
      int *Color1  = new int[npart+1];
      int *Color2  = new int[npart+1];
      double *px = new double[npart+1];
      double *py = new double[npart+1];
      double *pz = new double[npart+1];
      double *E = new double[npart+1];
      double *m = new double[npart+1];
      double *lifetime = new double[npart+1];
      double *spin = new double[npart+1];
      
      TLorentzVector **v = new TLorentzVector*[npart+1];
      TLorentzVector v_muon1, v_muon2, v_muon3, v_muon4, v_muon5, v_mu_nu1, v_mu_nu2, v_mu_nu3, v_mu_nu4, v_mu_nu5;
      TLorentzVector v_elec1, v_elec2, v_elec3, v_elec4, v_elec5, v_elec_nu1, v_elec_nu2, v_elec_nu3, v_elec_nu4, v_elec_nu5;
      TLorentzVector v_w_boson1, v_w_boson2, v_w_boson3, v_w_boson4, v_w_boson5;
      TLorentzVector v_bq1, v_bq2, v_bq3, v_bq4;
      TLorentzVector v_top1, v_top2, v_top3, v_top4, v_top5;
      TLorentzVector v_top_lep, v_top_alep, v_lep_nu, v_top_bJet, v_tbar_lep, v_top_jet, v_tbar_jet, v_tbar_nu, v_top_nu, v_wminus_lep;
      
      // in lhe first line is number 1, so fill unused array [0] with a crazy value;
      Id[0]= -99999;
      Status[0]= -99999;
      Mother1[0]= -99999;
      Mother2[0]= -99999;
      Color1[0]= -99999;
      Color2[0]= -99999;
      px[0]= -99999;
      py[0]= -99999;
      pz[0]= -99999;
      E[0]= -99999;
      m[0]= -99999;
      lifetime[0]= -99999;
      spin[0]= -99999;
     for (int i=1 ; i<npart+1 ; i++) { //start at one
        ff >> Id[i] >> Status[i] >> Mother1[i] >> Mother2[i] >> Color1[i] >> Color2[i]
           >> px[i] >> py[i] >> pz[i] >> E[i] >> m[i] >> lifetime[i] >> spin[i] ;
        buffer << Id[i] << " " << Status[i] << " " << std::endl;
        line++;
        v[i] = new TLorentzVector(px[i], py[i], pz[i], E[i]);
        if (Status[i]==-1) continue; // status -1 = initial quark ==> skip
        if (Id[i]==6)  top=i;
        if (Id[i]==-6) topbar=i;
        if (Id[i]>6000000) zprime=i;


//------------------top and tbar lep, bJets-----------------------------------
          
          if (Id[Mother1[i]]==24) 
          { // mother = W+ coming from top
          if ( Id[i]==-11 || Id[i]==-13 ) 
          { // charged anti leptons
          v_top_alep.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
          hPt_alep->Fill(v_top_alep.Pt());
          hEta_alep->Fill(v_top_alep.Eta());
          hPhi_alep->Fill(v_top_alep.Phi());
          hPx_alep->Fill(v_top_alep.Px());
          hPy_alep->Fill(v_top_alep.Py());
          n_alep++;
          }
          if ( Id[i]==12 || Id[i]==14 )
          {//neutrino
          v_lep_nu.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
          n_lep_nu++;
          }
          }
          
          if ( Id[Mother1[i]]==6 ) 
          { // mother = top
          if (Id[i]==5) 
          { // bjets
          v_top_bJet.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
          n_topbjets++;
          }
          }
          
          
          if ( Id[Mother1[i]==-24]) 
          {
          if ( Id[i]==11 || Id[i]==13 ) 
          { // charged leptons
          v_tbar_lep.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
          hPt_lep->Fill(v_tbar_lep.Pt());
          hEta_lep->Fill(v_tbar_lep.Eta());
          hPhi_lep->Fill(v_tbar_lep.Phi());
          hPx_lep->Fill(v_tbar_lep.Px());
          hPy_lep->Fill(v_tbar_lep.Py());
          n_lep++;
          }
          if ( Id[i]==-12 || Id[i]==-14 ) 
          {// anit neutrino
          v_tbar_nu.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
          n_tbar_nu++;
          }  
          }
          
          if ( abs(Id[Mother1[i]])== 6 ) 
          { // mother = tbar
          if ( abs(Id[i]) == 5) 
          { // bbarjets
          v_top_bJet.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
          hPt_bJet->Fill(v_top_bJet.Pt());
          hEta_bJet->Fill(v_top_bJet.Eta());
          hPhi_bJet->Fill(v_top_bJet.Phi());
          n_bJets++;
          }
          }

//-----------------------------------W boson reconstruction----------------------------------//

	if (abs(Id[Mother1[i]])==6)
	{ //Mother top
	  if (abs(Id[i])==5) { // bjets
	         if (bq[0] == -1) bq[0]=i;
            else if (bq[1]==-1) bq[1]=i;
            else if (bq[2]==-1) bq[2]=i;
            else if (bq[3]==-1) bq[3]=i;
            else cout << "ERROR : more than 4 b quarks" << endl;
           }
	}
	
	if ( abs(Id[Mother1[i]])==24)
	{
          if ( abs(Id[i])==11 ) 
          { // charged leptons
          	if( elec[0]== -1 ) elec[0]=i;
          	else if( elec[1]== -1 ) elec[1]=i;
		else if( elec[2]== -1 ) elec[2]=i;
		else if( elec[3]== -1 ) elec[3]=i;
		else if( elec[4]== -1 ) elec[4]=i;
		else cout<<"more than 5 electrons"<<endl;
		n_elec++;
	  }
	  if ( abs(Id[i])==12 ) 
          { // charged leptons neutrinos
          	if( elec_nu[0]== -1 ) elec_nu[0]=i;
          	else if( elec_nu[1]== -1 ) elec_nu[1]=i;
		else if( elec_nu[2]== -1 ) elec_nu[2]=i;
		else if( elec_nu[3]== -1 ) elec_nu[3]=i;
		else if( elec_nu[4]== -1 ) elec_nu[4]=i;
		else cout<<"more than 5 electron neutrinos"<<endl;
		n_elec_nu++;
	  }
	  if ( abs(Id[i])==13 ) 
          { // charged leptons
          	if( muon[0]== -1 ) muon[0]=i;
          	else if( muon[1]== -1 ) muon[1]=i;
		else if( muon[2]== -1 ) muon[2]=i;
		else if( muon[3]== -1 ) muon[3]=i;
		else if( muon[4]== -1 ) muon[4]=i;
		else cout<<"more than 5 muons"<<endl;
		n_muon++;
	  }
	  if ( abs(Id[i])==14 ) 
          { // charged leptons neutrinos
          	if( muon_nu[0]== -1 ) muon_nu[0]=i;
          	else if( muon_nu[1]== -1 ) muon_nu[1]=i;
		else if( muon_nu[2]== -1 ) muon_nu[2]=i;
		else if( muon_nu[3]== -1 ) muon_nu[3]=i;
		else if( muon_nu[4]== -1 ) muon_nu[4]=i;
		else cout<<"more than 5 muon neutrinos"<<endl;
		n_muon_nu++;
	  }
	}
} //loop of i

if (n_muon==5)
{
double muon_pt[5] = {v[muon[0]]->Pt(),v[muon[1]]->Pt(),v[muon[2]]->Pt(),v[muon[3]]->Pt(),v[muon[4]]->Pt()};	
        std::sort(muon_pt,muon_pt+5);
double muon_eta[5] = {v[muon[0]]->Eta(),v[muon[1]]->Eta(),v[muon[2]]->Eta(),v[muon[3]]->Eta(),v[muon[3]]->Eta()};	
        std::sort(muon_eta,muon_eta+5);
double muon_phi[5] = {v[muon[0]]->Phi(),v[muon[1]]->Phi(),v[muon[2]]->Phi(),v[muon[3]]->Phi(),v[muon[3]]->Phi()};	
        std::sort(muon_phi,muon_phi+5);
        
double bq_pt[4] = {v[bq[0]]->Pt(),v[bq[1]]->Pt(),v[bq[2]]->Pt(),v[bq[3]]->Pt()};	
        std::sort(bq_pt,bq_pt+4);
double bq_eta[4] = {v[bq[0]]->Eta(),v[bq[1]]->Eta(),v[bq[2]]->Eta(),v[bq[3]]->Eta()};	
        std::sort(bq_eta,bq_eta+4);
double bq_phi[4] = {v[bq[0]]->Phi(),v[bq[1]]->Phi(),v[bq[2]]->Phi(),v[bq[3]]->Phi()};	
        std::sort(bq_phi,bq_phi+4);
        
      v_bq1.SetPxPyPzE(v[bq[0]]->Px(),v[bq[0]]->Py(),v[bq[0]]->Pz(),v[bq[0]]->E());
      v_bq2.SetPxPyPzE(v[bq[1]]->Px(),v[bq[1]]->Py(),v[bq[1]]->Pz(),v[bq[1]]->E());
      v_bq3.SetPxPyPzE(v[bq[2]]->Px(),v[bq[2]]->Py(),v[bq[2]]->Pz(),v[bq[2]]->E());
      v_bq4.SetPxPyPzE(v[bq[3]]->Px(),v[bq[3]]->Py(),v[bq[3]]->Pz(),v[bq[3]]->E());
                
      v_muon1.SetPxPyPzE(v[muon[0]]->Px(),v[muon[0]]->Py(),v[muon[0]]->Pz(),v[muon[0]]->E());
      v_muon2.SetPxPyPzE(v[muon[1]]->Px(),v[muon[1]]->Py(),v[muon[1]]->Pz(),v[muon[1]]->E());
      v_muon3.SetPxPyPzE(v[muon[2]]->Px(),v[muon[2]]->Py(),v[muon[2]]->Pz(),v[muon[2]]->E());
      v_muon4.SetPxPyPzE(v[muon[3]]->Px(),v[muon[3]]->Py(),v[muon[3]]->Pz(),v[muon[3]]->E());
      v_muon5.SetPxPyPzE(v[muon[4]]->Px(),v[muon[4]]->Py(),v[muon[4]]->Pz(),v[muon[4]]->E());
      
      v_mu_nu1.SetPxPyPzE(v[muon_nu[0]]->Px(),v[muon_nu[0]]->Py(),v[muon_nu[0]]->Pz(),v[muon_nu[0]]->E());
      v_mu_nu2.SetPxPyPzE(v[muon_nu[1]]->Px(),v[muon_nu[1]]->Py(),v[muon_nu[1]]->Pz(),v[muon_nu[1]]->E());
      v_mu_nu3.SetPxPyPzE(v[muon_nu[2]]->Px(),v[muon_nu[2]]->Py(),v[muon_nu[2]]->Pz(),v[muon_nu[2]]->E());
      v_mu_nu4.SetPxPyPzE(v[muon_nu[3]]->Px(),v[muon_nu[3]]->Py(),v[muon_nu[3]]->Pz(),v[muon_nu[3]]->E());
      v_mu_nu5.SetPxPyPzE(v[muon_nu[4]]->Px(),v[muon_nu[4]]->Py(),v[muon_nu[4]]->Pz(),v[muon_nu[4]]->E());
      
      cout<<"bq_pt:"<<bq_pt[0]<<", "<<bq_pt[1]<<", "<<bq_pt[2]<<", "<<bq_pt[3]<<endl;
      cout<<"muon_pt :  "<<muon_pt[0]<<",  "<<muon_pt[1]<<",  "<<muon_pt[2]<<",  "<<muon_pt[3]<<",  "<<muon_pt[4]<<endl;
            
        // w boson kinematics:
        v_w_boson1 = v_muon1+v_mu_nu1;
        v_w_boson2 = v_muon2+v_mu_nu2;
        v_w_boson3 = v_muon3+v_mu_nu3;
        v_w_boson4 = v_muon4+v_mu_nu4;
        v_w_boson5 = v_muon5+v_mu_nu5;
        
        hPt_wboson->Fill(v_w_boson1.Pt());
        hEta_wboson->Fill(v_w_boson1.Eta());
        hPhi_wboson->Fill(v_w_boson1.Phi());
        
        hM_wboson1->Fill(v_w_boson1.M());
        hM_wboson2->Fill(v_w_boson2.M());
        hM_wboson3->Fill(v_w_boson3.M());
        hM_wboson4->Fill(v_w_boson4.M());
        hM_wboson5->Fill(v_w_boson5.M());
        
        // top quark kinematics
        v_top1= v_w_boson1+v_bq1;
        v_top2= v_w_boson2+v_bq2;
        v_top3= v_w_boson3+v_bq3;
        v_top4= v_w_boson4+v_bq4;
        //v_top5= v_w_boson5+v_bq4;
        
        hM_top1->Fill(v_top1.M());
        hM_top2->Fill(v_top2.M());
        hM_top3->Fill(v_top3.M());
        hM_top4->Fill(v_top4.M());
        //hM_top5->Fill(v_top5.M());
}

if (n_muon==4 )
{
double muon_pt[5] = {v[muon[0]]->Pt(),v[muon[1]]->Pt(),v[muon[2]]->Pt(),v[muon[3]]->Pt()};	
        std::sort(muon_pt,muon_pt+4);
double muon_eta[5] = {v[muon[0]]->Eta(),v[muon[1]]->Eta(),v[muon[2]]->Eta(),v[muon[3]]->Eta()};	
        std::sort(muon_eta,muon_eta+4);
double muon_phi[5] = {v[muon[0]]->Phi(),v[muon[1]]->Phi(),v[muon[2]]->Phi(),v[muon[3]]->Phi()};	
        std::sort(muon_phi,muon_phi+4);

double elec_pt[5] = {v[elec[0]]->Pt()};	
        std::sort(elec_pt,elec_pt+1);
double elec_eta[5] = {v[elec[0]]->Eta()};	
        std::sort(elec_eta,elec_eta+1);
double elec_phi[5] = {v[elec[0]]->Phi()};	
        std::sort(elec_phi,elec_phi+1);
        
double bq_pt[4] = {v[bq[0]]->Pt(),v[bq[1]]->Pt(),v[bq[2]]->Pt(),v[bq[3]]->Pt()};	
        std::sort(bq_pt,bq_pt+4);
double bq_eta[4] = {v[bq[0]]->Eta(),v[bq[1]]->Eta(),v[bq[2]]->Eta(),v[bq[3]]->Eta()};	
        std::sort(bq_eta,bq_eta+4);
double bq_phi[4] = {v[bq[0]]->Phi(),v[bq[1]]->Phi(),v[bq[2]]->Phi(),v[bq[3]]->Phi()};	
        std::sort(bq_phi,bq_phi+4);
        
      v_bq1.SetPxPyPzE(v[bq[0]]->Px(),v[bq[0]]->Py(),v[bq[0]]->Pz(),v[bq[0]]->E());
      v_bq2.SetPxPyPzE(v[bq[1]]->Px(),v[bq[1]]->Py(),v[bq[1]]->Pz(),v[bq[1]]->E());
      v_bq3.SetPxPyPzE(v[bq[2]]->Px(),v[bq[2]]->Py(),v[bq[2]]->Pz(),v[bq[2]]->E());
      v_bq4.SetPxPyPzE(v[bq[3]]->Px(),v[bq[3]]->Py(),v[bq[3]]->Pz(),v[bq[3]]->E());
                
      v_muon1.SetPxPyPzE(v[muon[0]]->Px(),v[muon[0]]->Py(),v[muon[0]]->Pz(),v[muon[0]]->E());
      v_muon2.SetPxPyPzE(v[muon[1]]->Px(),v[muon[1]]->Py(),v[muon[1]]->Pz(),v[muon[1]]->E());
      v_muon3.SetPxPyPzE(v[muon[2]]->Px(),v[muon[2]]->Py(),v[muon[2]]->Pz(),v[muon[2]]->E());
      v_muon4.SetPxPyPzE(v[muon[3]]->Px(),v[muon[3]]->Py(),v[muon[3]]->Pz(),v[muon[3]]->E());
      
      v_mu_nu1.SetPxPyPzE(v[muon_nu[0]]->Px(),v[muon_nu[0]]->Py(),v[muon_nu[0]]->Pz(),v[muon_nu[0]]->E());
      v_mu_nu2.SetPxPyPzE(v[muon_nu[1]]->Px(),v[muon_nu[1]]->Py(),v[muon_nu[1]]->Pz(),v[muon_nu[1]]->E());
      v_mu_nu3.SetPxPyPzE(v[muon_nu[2]]->Px(),v[muon_nu[2]]->Py(),v[muon_nu[2]]->Pz(),v[muon_nu[2]]->E());
      v_mu_nu4.SetPxPyPzE(v[muon_nu[3]]->Px(),v[muon_nu[3]]->Py(),v[muon_nu[3]]->Pz(),v[muon_nu[3]]->E());
      
      v_elec1.SetPxPyPzE(v[elec[0]]->Px(),v[elec[0]]->Py(),v[elec[0]]->Pz(),v[elec[0]]->E());
      v_elec_nu1.SetPxPyPzE(v[elec_nu[0]]->Px(),v[elec_nu[0]]->Py(),v[elec_nu[0]]->Pz(),v[elec_nu[0]]->E());
      
      cout<<"bq_pt:"<<bq_pt[0]<<", "<<bq_pt[1]<<", "<<bq_pt[2]<<", "<<bq_pt[3]<<endl;
      cout<<"muon_pt :  "<<muon_pt[0]<<",  "<<muon_pt[1]<<",  "<<muon_pt[2]<<",  "<<muon_pt[3]<<endl;
      cout<<"elec_pt :  "<<elec_pt[0]<<", "<<endl;
            
        // w boson kinematics:
        v_w_boson1 = v_muon1+v_mu_nu1;
        v_w_boson2 = v_muon2+v_mu_nu2;
        v_w_boson3 = v_muon3+v_mu_nu3;
        v_w_boson4 = v_muon4+v_mu_nu4;
        v_w_boson5 = v_elec1+v_elec_nu1;
        
        hPt_wboson->Fill(v_w_boson1.Pt());
        hEta_wboson->Fill(v_w_boson1.Eta());
        hPhi_wboson->Fill(v_w_boson1.Phi());
        
        hPt_wboson->Fill(v_w_boson1.Pt());
        hM_wboson2->Fill(v_w_boson2.M());
        hM_wboson3->Fill(v_w_boson3.M());
        hM_wboson4->Fill(v_w_boson4.M());
        hM_wboson5->Fill(v_w_boson5.M());
        
        // top quark kinematics
        v_top1= v_w_boson1+v_bq1;
        v_top2= v_w_boson2+v_bq2;
        v_top3= v_w_boson3+v_bq3;
        v_top4= v_w_boson4+v_bq4;
        //v_top5= v_w_boson5+v_bq2; 
        
        hM_top1->Fill(v_top1.M());
        hM_top2->Fill(v_top2.M());
        hM_top3->Fill(v_top3.M());
        hM_top4->Fill(v_top4.M());
        //hM_top5->Fill(v_top5.M());
}

if (n_muon==3 )
{
double muon_pt[5] = {v[muon[0]]->Pt(),v[muon[1]]->Pt(),v[muon[2]]->Pt()};	
        std::sort(muon_pt,muon_pt+3);
double muon_eta[5] = {v[muon[0]]->Eta(),v[muon[1]]->Eta(),v[muon[2]]->Eta()};	
        std::sort(muon_eta,muon_eta+3);
double muon_phi[5] = {v[muon[0]]->Phi(),v[muon[1]]->Phi(),v[muon[2]]->Phi()};	
        std::sort(muon_phi,muon_phi+3);
        
double elec_pt[5] = {v[elec[0]]->Pt(), v[elec[1]]->Pt()};	
        std::sort(elec_pt,elec_pt+2);
double elec_eta[5] = {v[elec[0]]->Eta(), v[elec[1]]->Eta()};	
        std::sort(elec_eta,elec_eta+2);
double elec_phi[5] = {v[elec[0]]->Phi(), v[elec[1]]->Phi()};	
        std::sort(elec_phi,elec_phi+2);
        
double bq_pt[4] = {v[bq[0]]->Pt(),v[bq[1]]->Pt(),v[bq[2]]->Pt(),v[bq[3]]->Pt()};	
        std::sort(bq_pt,bq_pt+4);
double bq_eta[4] = {v[bq[0]]->Eta(),v[bq[1]]->Eta(),v[bq[2]]->Eta(),v[bq[3]]->Eta()};	
        std::sort(bq_eta,bq_eta+4);
double bq_phi[4] = {v[bq[0]]->Phi(),v[bq[1]]->Phi(),v[bq[2]]->Phi(),v[bq[3]]->Phi()};	
        std::sort(bq_phi,bq_phi+4);
        
      v_bq1.SetPxPyPzE(v[bq[0]]->Px(),v[bq[0]]->Py(),v[bq[0]]->Pz(),v[bq[0]]->E());
      v_bq2.SetPxPyPzE(v[bq[1]]->Px(),v[bq[1]]->Py(),v[bq[1]]->Pz(),v[bq[1]]->E());
      v_bq3.SetPxPyPzE(v[bq[2]]->Px(),v[bq[2]]->Py(),v[bq[2]]->Pz(),v[bq[2]]->E());
      v_bq4.SetPxPyPzE(v[bq[3]]->Px(),v[bq[3]]->Py(),v[bq[3]]->Pz(),v[bq[3]]->E());
        
      v_muon1.SetPxPyPzE(v[muon[0]]->Px(),v[muon[0]]->Py(),v[muon[0]]->Pz(),v[muon[0]]->E());
      v_muon2.SetPxPyPzE(v[muon[1]]->Px(),v[muon[1]]->Py(),v[muon[1]]->Pz(),v[muon[1]]->E());
      v_muon3.SetPxPyPzE(v[muon[2]]->Px(),v[muon[2]]->Py(),v[muon[2]]->Pz(),v[muon[2]]->E());
      
      v_mu_nu1.SetPxPyPzE(v[muon_nu[0]]->Px(),v[muon_nu[0]]->Py(),v[muon_nu[0]]->Pz(),v[muon_nu[0]]->E());
      v_mu_nu2.SetPxPyPzE(v[muon_nu[1]]->Px(),v[muon_nu[1]]->Py(),v[muon_nu[1]]->Pz(),v[muon_nu[1]]->E());
      v_mu_nu3.SetPxPyPzE(v[muon_nu[2]]->Px(),v[muon_nu[2]]->Py(),v[muon_nu[2]]->Pz(),v[muon_nu[2]]->E());
      
      v_elec1.SetPxPyPzE(v[elec[0]]->Px(),v[elec[0]]->Py(),v[elec[0]]->Pz(),v[elec[0]]->E());
      v_elec2.SetPxPyPzE(v[elec[1]]->Px(),v[elec[1]]->Py(),v[elec[1]]->Pz(),v[elec[1]]->E());
      v_elec_nu1.SetPxPyPzE(v[elec_nu[0]]->Px(),v[elec_nu[0]]->Py(),v[elec_nu[0]]->Pz(),v[elec_nu[0]]->E());
      v_elec_nu2.SetPxPyPzE(v[elec_nu[1]]->Px(),v[elec_nu[1]]->Py(),v[elec_nu[1]]->Pz(),v[elec_nu[1]]->E());
      
      cout<<"bq_pt:"<<bq_pt[0]<<", "<<bq_pt[1]<<", "<<bq_pt[2]<<", "<<bq_pt[3]<<endl;
      cout<<"muon_pt :  "<<muon_pt[0]<<",  "<<muon_pt[1]<<",  "<<muon_pt[2]<<endl;
      cout<<"elec_pt :  "<<elec_pt[0]<<", "<<elec_pt[1]<<", "<<endl;
            
        // w boson kinematics:
        v_w_boson1 = v_muon1+v_mu_nu1;
        v_w_boson2 = v_muon2+v_mu_nu2;
        v_w_boson3 = v_muon3+v_mu_nu3;
        v_w_boson4 = v_elec1+v_elec_nu1;
        v_w_boson5 = v_elec2+v_elec_nu2;
        
        hPt_wboson->Fill(v_w_boson1.Pt());
        hEta_wboson->Fill(v_w_boson1.Eta());
        hPhi_wboson->Fill(v_w_boson1.Phi());
        
        hPt_wboson->Fill(v_w_boson1.Pt());
        hM_wboson2->Fill(v_w_boson2.M());
        hM_wboson3->Fill(v_w_boson3.M());
        hM_wboson4->Fill(v_w_boson4.M());
        hM_wboson5->Fill(v_w_boson5.M());
        
        // top quark kinematics
        v_top1= v_w_boson1+v_bq1;
        v_top2= v_w_boson2+v_bq2;
        v_top3= v_w_boson3+v_bq3;
        v_top4= v_w_boson4+v_bq4;
        //v_top5= v_w_boson5+v_bq3;
        
        hM_top1->Fill(v_top1.M());
        hM_top2->Fill(v_top2.M());
        hM_top3->Fill(v_top3.M());
        hM_top4->Fill(v_top4.M());
        //hM_top5->Fill(v_top5.M());
}

if (n_muon==2 )
{
double muon_pt[5] = {v[muon[0]]->Pt(),v[muon[1]]->Pt()};	
        std::sort(muon_pt,muon_pt+2);
double muon_eta[5] = {v[muon[0]]->Eta(),v[muon[1]]->Eta()};	
        std::sort(muon_eta,muon_eta+2);
double muon_phi[5] = {v[muon[0]]->Phi(),v[muon[1]]->Phi()};	
        std::sort(muon_phi,muon_phi+2);

double elec_pt[5] = {v[elec[0]]->Pt(), v[elec[1]]->Pt(), v[elec[2]]->Pt()};	
        std::sort(elec_pt,elec_pt+3);
double elec_eta[5] = {v[elec[0]]->Eta(), v[elec[1]]->Eta(), v[elec[2]]->Eta()};	
        std::sort(elec_eta,elec_eta+3);
double elec_phi[5] = {v[elec[0]]->Phi(), v[elec[1]]->Phi(), v[elec[2]]->Phi()};	
        std::sort(elec_phi,elec_phi+3);
        
double bq_pt[4] = {v[bq[0]]->Pt(),v[bq[1]]->Pt(),v[bq[2]]->Pt(),v[bq[3]]->Pt()};	
        std::sort(bq_pt,bq_pt+4);
double bq_eta[4] = {v[bq[0]]->Eta(),v[bq[1]]->Eta(),v[bq[2]]->Eta(),v[bq[3]]->Eta()};	
        std::sort(bq_eta,bq_eta+4);
double bq_phi[4] = {v[bq[0]]->Phi(),v[bq[1]]->Phi(),v[bq[2]]->Phi(),v[bq[3]]->Phi()};	
        std::sort(bq_phi,bq_phi+4);
        
      v_bq1.SetPxPyPzE(v[bq[0]]->Px(),v[bq[0]]->Py(),v[bq[0]]->Pz(),v[bq[0]]->E());
      v_bq2.SetPxPyPzE(v[bq[1]]->Px(),v[bq[1]]->Py(),v[bq[1]]->Pz(),v[bq[1]]->E());
      v_bq3.SetPxPyPzE(v[bq[2]]->Px(),v[bq[2]]->Py(),v[bq[2]]->Pz(),v[bq[2]]->E());
      v_bq4.SetPxPyPzE(v[bq[3]]->Px(),v[bq[3]]->Py(),v[bq[3]]->Pz(),v[bq[3]]->E());
        
      v_muon1.SetPxPyPzE(v[muon[0]]->Px(),v[muon[0]]->Py(),v[muon[0]]->Pz(),v[muon[0]]->E());
      v_muon2.SetPxPyPzE(v[muon[1]]->Px(),v[muon[1]]->Py(),v[muon[1]]->Pz(),v[muon[1]]->E());
      v_mu_nu1.SetPxPyPzE(v[muon_nu[0]]->Px(),v[muon_nu[0]]->Py(),v[muon_nu[0]]->Pz(),v[muon_nu[0]]->E());
      v_mu_nu2.SetPxPyPzE(v[muon_nu[1]]->Px(),v[muon_nu[1]]->Py(),v[muon_nu[1]]->Pz(),v[muon_nu[1]]->E());
      
      v_elec1.SetPxPyPzE(v[elec[0]]->Px(),v[elec[0]]->Py(),v[elec[0]]->Pz(),v[elec[0]]->E());
      v_elec2.SetPxPyPzE(v[elec[1]]->Px(),v[elec[1]]->Py(),v[elec[1]]->Pz(),v[elec[1]]->E());
      v_elec3.SetPxPyPzE(v[elec[2]]->Px(),v[elec[2]]->Py(),v[elec[2]]->Pz(),v[elec[2]]->E());
      v_elec_nu1.SetPxPyPzE(v[elec_nu[0]]->Px(),v[elec_nu[0]]->Py(),v[elec_nu[0]]->Pz(),v[elec_nu[0]]->E());
      v_elec_nu2.SetPxPyPzE(v[elec_nu[1]]->Px(),v[elec_nu[1]]->Py(),v[elec_nu[1]]->Pz(),v[elec_nu[1]]->E());
      v_elec_nu3.SetPxPyPzE(v[elec_nu[2]]->Px(),v[elec_nu[2]]->Py(),v[elec_nu[2]]->Pz(),v[elec_nu[2]]->E());
      
      cout<<"bq_pt:"<<bq_pt[0]<<", "<<bq_pt[1]<<", "<<bq_pt[2]<<", "<<bq_pt[3]<<endl;
      cout<<"muon_pt :  "<<muon_pt[0]<<",  "<<muon_pt[1]<<endl;
      cout<<"elec_pt :  "<<elec_pt[0]<<", "<<elec_pt[1]<<", "<<elec_pt[2]<<", "<<endl;
            
        // w boson kinematics:
        v_w_boson1 = v_muon1+v_mu_nu1;
        v_w_boson2 = v_muon2+v_mu_nu2;
        v_w_boson3 = v_elec1+v_elec_nu1;
        v_w_boson4 = v_elec2+v_elec_nu2;
        v_w_boson5 = v_elec3+v_elec_nu3;
        
        hPt_wboson->Fill(v_w_boson1.Pt());
        hEta_wboson->Fill(v_w_boson1.Eta());
        hPhi_wboson->Fill(v_w_boson1.Phi());
        
        hPt_wboson->Fill(v_w_boson1.Pt());
        hM_wboson2->Fill(v_w_boson2.M());
        hM_wboson3->Fill(v_w_boson3.M());
        hM_wboson4->Fill(v_w_boson4.M());
        hM_wboson5->Fill(v_w_boson5.M());
        
        // top quark kinematics
        v_top1= v_w_boson1+v_bq1;
        v_top2= v_w_boson2+v_bq2;
        v_top3= v_w_boson3+v_bq3;
        v_top4= v_w_boson4+v_bq4;
        
        hM_top1->Fill(v_top1.M());
        hM_top2->Fill(v_top2.M());
        hM_top3->Fill(v_top3.M());
        hM_top4->Fill(v_top4.M());
}

if (n_muon==1 )
{
double muon_pt[5] = {v[muon[0]]->Pt()};	
        std::sort(muon_pt,muon_pt+1);
double muon_eta[5] = {v[muon[0]]->Eta()};	
        std::sort(muon_eta,muon_eta+1);
double muon_phi[5] = {v[muon[0]]->Phi()};	
        std::sort(muon_phi,muon_phi+1);

double elec_pt[5] = {v[elec[0]]->Pt(), v[elec[1]]->Pt(), v[elec[2]]->Pt(), v[elec[3]]->Pt()};	
        std::sort(elec_pt,elec_pt+4);
double elec_eta[5] = {v[elec[0]]->Eta(), v[elec[1]]->Eta(), v[elec[2]]->Eta(), v[elec[3]]->Eta()};	
        std::sort(elec_eta,elec_eta+4);
double elec_phi[5] = {v[elec[0]]->Phi(), v[elec[1]]->Phi(), v[elec[2]]->Phi(), v[elec[3]]->Phi()};	
        std::sort(elec_phi,elec_phi+4);
        
double bq_pt[4] = {v[bq[0]]->Pt(),v[bq[1]]->Pt(),v[bq[2]]->Pt(),v[bq[3]]->Pt()};	
        std::sort(bq_pt,bq_pt+4);
double bq_eta[4] = {v[bq[0]]->Eta(),v[bq[1]]->Eta(),v[bq[2]]->Eta(),v[bq[3]]->Eta()};	
        std::sort(bq_eta,bq_eta+4);
double bq_phi[4] = {v[bq[0]]->Phi(),v[bq[1]]->Phi(),v[bq[2]]->Phi(),v[bq[3]]->Phi()};	
        std::sort(bq_phi,bq_phi+4);
        
      v_bq1.SetPxPyPzE(v[bq[0]]->Px(),v[bq[0]]->Py(),v[bq[0]]->Pz(),v[bq[0]]->E());
      v_bq2.SetPxPyPzE(v[bq[1]]->Px(),v[bq[1]]->Py(),v[bq[1]]->Pz(),v[bq[1]]->E());
      v_bq3.SetPxPyPzE(v[bq[2]]->Px(),v[bq[2]]->Py(),v[bq[2]]->Pz(),v[bq[2]]->E());
      v_bq4.SetPxPyPzE(v[bq[3]]->Px(),v[bq[3]]->Py(),v[bq[3]]->Pz(),v[bq[3]]->E());
        
      v_muon1.SetPxPyPzE(v[muon[0]]->Px(),v[muon[0]]->Py(),v[muon[0]]->Pz(),v[muon[0]]->E());
      v_mu_nu1.SetPxPyPzE(v[muon_nu[0]]->Px(),v[muon_nu[0]]->Py(),v[muon_nu[0]]->Pz(),v[muon_nu[0]]->E());
      
      v_elec1.SetPxPyPzE(v[elec[0]]->Px(),v[elec[0]]->Py(),v[elec[0]]->Pz(),v[elec[0]]->E());
      v_elec2.SetPxPyPzE(v[elec[1]]->Px(),v[elec[1]]->Py(),v[elec[1]]->Pz(),v[elec[1]]->E());
      v_elec3.SetPxPyPzE(v[elec[2]]->Px(),v[elec[2]]->Py(),v[elec[2]]->Pz(),v[elec[2]]->E());
      v_elec4.SetPxPyPzE(v[elec[3]]->Px(),v[elec[3]]->Py(),v[elec[3]]->Pz(),v[elec[3]]->E());
      v_elec_nu1.SetPxPyPzE(v[elec_nu[0]]->Px(),v[elec_nu[0]]->Py(),v[elec_nu[0]]->Pz(),v[elec_nu[0]]->E());
      v_elec_nu2.SetPxPyPzE(v[elec_nu[1]]->Px(),v[elec_nu[1]]->Py(),v[elec_nu[1]]->Pz(),v[elec_nu[1]]->E());
      v_elec_nu3.SetPxPyPzE(v[elec_nu[2]]->Px(),v[elec_nu[2]]->Py(),v[elec_nu[2]]->Pz(),v[elec_nu[2]]->E());
      v_elec_nu4.SetPxPyPzE(v[elec_nu[3]]->Px(),v[elec_nu[3]]->Py(),v[elec_nu[3]]->Pz(),v[elec_nu[3]]->E());
      
      cout<<"bq_pt:"<<bq_pt[0]<<", "<<bq_pt[1]<<", "<<bq_pt[2]<<", "<<bq_pt[3]<<endl;
      cout<<"muon_pt :  "<<muon_pt[0]<<endl;
      cout<<"elec_pt :  "<<elec_pt[0]<<", "<<elec_pt[1]<<", "<<elec_pt[2]<<", "<<elec_pt[3]<<", "<<endl;
            
        // w boson kinematics:
        v_w_boson1 = v_muon1+v_mu_nu1;
        v_w_boson2 = v_elec1+v_elec_nu1;
        v_w_boson3 = v_elec2+v_elec_nu2;
        v_w_boson4 = v_elec3+v_elec_nu3;
        v_w_boson5 = v_elec4+v_elec_nu4;
        
        hPt_wboson->Fill(v_w_boson1.Pt());
        hEta_wboson->Fill(v_w_boson1.Eta());
        hPhi_wboson->Fill(v_w_boson1.Phi());
        
        hPt_wboson->Fill(v_w_boson1.Pt());
        hM_wboson2->Fill(v_w_boson2.M());
        hM_wboson3->Fill(v_w_boson3.M());
        hM_wboson4->Fill(v_w_boson4.M());
        hM_wboson5->Fill(v_w_boson5.M());
        
        // top quark kinematics:
        v_top1= v_w_boson1+v_bq1;
        v_top2= v_w_boson2+v_bq2;
        v_top3= v_w_boson3+v_bq3;
        v_top4= v_w_boson4+v_bq4;
        
        hM_top1->Fill(v_top1.M());
        hM_top2->Fill(v_top2.M());
        hM_top3->Fill(v_top3.M());
        hM_top4->Fill(v_top4.M());
}

if (n_elec==5 )
{
double elec_pt[5] = {v[elec[0]]->Pt(), v[elec[1]]->Pt(), v[elec[2]]->Pt(), v[elec[3]]->Pt(), v[elec[4]]->Pt()};	
        std::sort(elec_pt,elec_pt+5);
double elec_eta[5] = {v[elec[0]]->Eta(), v[elec[1]]->Eta(), v[elec[2]]->Eta(), v[elec[3]]->Eta(), v[elec[4]]->Eta()};	
        std::sort(elec_eta,elec_eta+5);
double elec_phi[5] = {v[elec[0]]->Phi(), v[elec[1]]->Phi(), v[elec[2]]->Phi(), v[elec[3]]->Phi(), v[elec[4]]->Phi()};	
        std::sort(elec_phi,elec_phi+5);
        
double bq_pt[4] = {v[bq[0]]->Pt(),v[bq[1]]->Pt(),v[bq[2]]->Pt(),v[bq[3]]->Pt()};	
        std::sort(bq_pt,bq_pt+4);
double bq_eta[4] = {v[bq[0]]->Eta(),v[bq[1]]->Eta(),v[bq[2]]->Eta(),v[bq[3]]->Eta()};	
        std::sort(bq_eta,bq_eta+4);
double bq_phi[4] = {v[bq[0]]->Phi(),v[bq[1]]->Phi(),v[bq[2]]->Phi(),v[bq[3]]->Phi()};	
        std::sort(bq_phi,bq_phi+4);
        
      v_bq1.SetPxPyPzE(v[bq[0]]->Px(),v[bq[0]]->Py(),v[bq[0]]->Pz(),v[bq[0]]->E());
      v_bq2.SetPxPyPzE(v[bq[1]]->Px(),v[bq[1]]->Py(),v[bq[1]]->Pz(),v[bq[1]]->E());
      v_bq3.SetPxPyPzE(v[bq[2]]->Px(),v[bq[2]]->Py(),v[bq[2]]->Pz(),v[bq[2]]->E());
      v_bq4.SetPxPyPzE(v[bq[3]]->Px(),v[bq[3]]->Py(),v[bq[3]]->Pz(),v[bq[3]]->E());
        
      v_elec1.SetPxPyPzE(v[elec[0]]->Px(),v[elec[0]]->Py(),v[elec[0]]->Pz(),v[elec[0]]->E());
      v_elec2.SetPxPyPzE(v[elec[1]]->Px(),v[elec[1]]->Py(),v[elec[1]]->Pz(),v[elec[1]]->E());
      v_elec3.SetPxPyPzE(v[elec[2]]->Px(),v[elec[2]]->Py(),v[elec[2]]->Pz(),v[elec[2]]->E());
      v_elec4.SetPxPyPzE(v[elec[3]]->Px(),v[elec[3]]->Py(),v[elec[3]]->Pz(),v[elec[3]]->E());
      v_elec5.SetPxPyPzE(v[elec[4]]->Px(),v[elec[4]]->Py(),v[elec[4]]->Pz(),v[elec[4]]->E());
      
      v_elec_nu1.SetPxPyPzE(v[elec_nu[0]]->Px(),v[elec_nu[0]]->Py(),v[elec_nu[0]]->Pz(),v[elec_nu[0]]->E());
      v_elec_nu2.SetPxPyPzE(v[elec_nu[1]]->Px(),v[elec_nu[1]]->Py(),v[elec_nu[1]]->Pz(),v[elec_nu[1]]->E());
      v_elec_nu3.SetPxPyPzE(v[elec_nu[2]]->Px(),v[elec_nu[2]]->Py(),v[elec_nu[2]]->Pz(),v[elec_nu[2]]->E());
      v_elec_nu4.SetPxPyPzE(v[elec_nu[3]]->Px(),v[elec_nu[3]]->Py(),v[elec_nu[3]]->Pz(),v[elec_nu[3]]->E());
      v_elec_nu5.SetPxPyPzE(v[elec_nu[4]]->Px(),v[elec_nu[4]]->Py(),v[elec_nu[4]]->Pz(),v[elec_nu[4]]->E());
            
      cout<<"bq_pt:"<<bq_pt[0]<<", "<<bq_pt[1]<<", "<<bq_pt[2]<<", "<<bq_pt[3]<<endl;
      cout<<"elec_pt :  "<<elec_pt[0]<<", "<<elec_pt[1]<<", "<<elec_pt[2]<<", "<<elec_pt[3]<<", "<<elec_pt[4]<<", "<<endl;
            
        // w boson kinematics:
        v_w_boson1 = v_elec1+v_elec_nu1;
        v_w_boson2 = v_elec2+v_elec_nu2;
        v_w_boson3 = v_elec3+v_elec_nu3;
        v_w_boson4 = v_elec4+v_elec_nu4;
        v_w_boson5 = v_elec5+v_elec_nu5;
        
        hPt_wboson->Fill(v_w_boson1.Pt());
        hEta_wboson->Fill(v_w_boson1.Eta());
        hPhi_wboson->Fill(v_w_boson1.Phi());
        
        hPt_wboson->Fill(v_w_boson1.Pt());
        hM_wboson2->Fill(v_w_boson2.M());
        hM_wboson3->Fill(v_w_boson3.M());
        hM_wboson4->Fill(v_w_boson4.M());
        hM_wboson5->Fill(v_w_boson5.M());
        
        // top quark kinematics
        v_top1= v_w_boson1+v_bq1;
        v_top2= v_w_boson2+v_bq2;
        v_top3= v_w_boson3+v_bq3;
        v_top4= v_w_boson4+v_bq4;
        
        hM_top1->Fill(v_top1.M());
        hM_top2->Fill(v_top2.M());
        hM_top3->Fill(v_top3.M());
        hM_top4->Fill(v_top4.M());
}

// --- end filling the histograms

      ff>>tt;
      line++;
      //if (event==100)  break;
      delete Id;
      delete Status;
      delete Mother1;
      delete Mother2;
      delete Color1;
      delete Color2;
      delete px;
      delete py;
      delete pz;
      delete E;
      delete m;
      delete lifetime;
      delete spin;
      for (int k=1 ; k<npart+1 ; delete v[k++]);

    }
}

  cout << " Total number of events --> " << event << endl;
  TFile *rootfile = new TFile(rootfname.c_str(),"recreate");
  
  hPt_lep->Write();
  hEta_lep->Write();
  hPhi_lep->Write();
  hPx_lep->Write();
  hPy_lep->Write();
  
  hPt_alep->Write();
  hEta_alep->Write();
  hPhi_alep->Write();
  hPx_alep->Write();
  hPy_alep->Write();
  
  hPt_bJet->Write();
  hEta_bJet->Write();
  hPhi_bJet->Write();
  
  hPt_wboson->Write();
  hEta_wboson->Write();
  hPhi_wboson->Write();
  
  hM_wboson1->Write();
  hM_wboson2->Write();
  hM_wboson3->Write();
  hM_wboson4->Write();
  hM_wboson5->Write();

  hM_top1->Write();
  hM_top2->Write();
  hM_top3->Write();
  hM_top4->Write();
  hM_top5->Write();

  rootfile->Close();

  cout << "Events with negative weight: " << negativeWeight << endl;
  cout << "lept decay = " << nlept*1.0/event << endl;
  cout << "hadr decay = " << nhadr*1.0/event << endl;
  cout << "semi decay = " << nsemi*1.0/event << endl;
  exit(0);

}
float DeltaPhi(float phi1, float phi2)
{
  float dphi = phi2 - phi1;
  if (fabs(dphi) > 3.14) dphi = 6.28 - fabs(dphi);
  return dphi;
}
float DeltaR(float eta1, float eta2, float phi1, float phi2)
{
  float deta = eta2 - eta1;
  float dphi = phi2 - phi1;
  if (fabs(dphi) > 3.14) dphi = 6.28 - fabs(dphi);
  float DELTAR = sqrt(pow(dphi,2)+pow(deta,2))*1.0;
  return DELTAR;
}

