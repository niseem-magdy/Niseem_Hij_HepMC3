//---------------------------------------------------------------
//---------------------------------------------------------------
//-------------------This Niseem Magdy Work 12:15:2023
//---------------------------------------------------------------
//---------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <map>
#include <string>
#include <math.h>

#include <random>
#include <chrono>

#include <csignal>
//---------------------------------------------------------------
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Print.h"
#include "HepMC3/Selector.h"
//--------------------------------------------------------------- 
#include "TROOT.h"
#include "TFile.h"

#include "TChain.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TMath.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TUnixSystem.h"
#include "TRandom3.h"

#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include "TClassTable.h"
#include "TSystem.h"
//--------------------------------------------------------------- 
#ifdef DEBUG
#define LogDebug(message) std::cout << "[DEBUG] " << message << std::endl
#else
#define LogDebug(message)
#endif
//---------------------------------------------------------------
using namespace HepMC3;
//---------------------------------------------------------------
const double FMtoMM = 1.0e-12;
const double Pi     = 3.141592653589793;
//---------------------------------------------------------------
extern "C" {
  void hijset_(float* efrm, char* frame, char* proj, char* targ, int* iap, int* izp, int* iat, int* izt);
}
#define HIJSET hijset_

extern "C" {
  void hijing_(char* FRAME, float* BMIN, float* BMAX);
}
#define HIJING hijing_

extern "C" {
  extern struct{ 
    int natt;
    int eatt;
    int jatt;
    int nt;
    int np;
    int n0;
    int n01;
    int n10;
    int n11;
  }himain1_;
}
#define himain1 himain1_

extern "C" {
  extern struct{ 
    int   katt[4][130000];
    float patt[4][130000];
    float vatt[4][130000];
  }himain2_;
}
#define himain2 himain2_

extern "C" {
  extern struct{
    float vatt[4][130000];
  }himain3_;
}
#define himain3 himain3_

extern "C" {
  extern struct{ 
    float  hipr1[100];
    int    ihpr2[50];
    float  hint1[100];
    int    ihnt2[50];
  }hiparnt_;
}
#define hiparnt hiparnt_

extern "C" {
  extern struct{
    int   nfp[15][300];
    float pp[15][300];
    int   nft[15][300];
    float pt[15][300];   
  }histrng_;
}
#define histrng histrng_

extern "C" {
  extern struct{
    float yp[300][3];
    float yt[300][3];
  }hijcrdn_;
}
#define hijcrdn hijcrdn_

extern "C" {
  extern struct{
    int nseed;
  }ranseed_;
}
#define ranseed ranseed_
//---------------------------------------------------------------
int main(int argc, char* argv[]){
  //-------------------------------------------------------------
  if (argc != 12) {
    std::cerr << "Usage: " << argv[0] << " Mode Decay FRAME N_EVENT IAP IZP IAT IZT EFRM BMIN BMAX" << std::endl;
    return 1;
  }
  //-------------------------------------------------------------
  int MODE    = std::stoi(argv[1]);//0-HepMC 1 for tree
  int DECAY   = std::stoi(argv[2]);//0 or 1
  int FRAMEID = std::stoi(argv[3]);//0-CMS
  int N_EVENT = std::stoi(argv[4]);
  int IAP     = std::stoi(argv[5]);
  int IZP     = std::stoi(argv[6]);
  int IAT     = std::stoi(argv[7]);
  int IZT     = std::stoi(argv[8]);
  float EFRM  = std::stof(argv[9]);
  float BMIN  = std::stof(argv[10]);
  float BMAX  = std::stof(argv[11]); 

  char FRAME[2][9];
  std::strcpy(FRAME[0], "CMS     ");
  std::strcpy(FRAME[1], "LAB     ");

  char PROJ[2][9];
  std::strcpy(PROJ[0], "P       ");
  std::strcpy(PROJ[1], "A       ");

  char TARG[2][9];
  std::strcpy(TARG[0],"P       ");
  std::strcpy(TARG[1],"A       ");
  
  int PROJID;
  if(IAP == 1){
    PROJID = 0;
  }else{
    PROJID = 1;
  }

  int TARGID;
  if(IAT == 1){
    TARGID = 0;
  }else{
    TARGID = 1;
  }

  float RIonP    = 1.2*pow(IAP,0.3333333333);
  float RIonT    = 1.2*pow(IAT,0.3333333333);
  float RPT      =  (RIonP+RIonT)/2.0;
  //--------------------------------------------------------------
  // Get the current time as a duration since the epoch
  auto currentTime = std::chrono::system_clock::now().time_since_epoch();
  auto currentTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime).count();
  // Seed the random number generator with the current time
  std::mt19937_64 generator(currentTimeMs);
  // Generate a random odd number
  std::uniform_int_distribution<int> distribution(1, 20000);
  int randomNumber = distribution(generator);
  if (randomNumber % 2 == 0) {
    randomNumber++; // Make sure it's odd
  }  
  std::cout << "Random odd number: " << randomNumber << std::endl;
  //--------------------------------------------------------------
  ranseed.nseed    = randomNumber;
  hiparnt.ihpr2[9] =0;
  hiparnt.ihpr2[20]=DECAY;
  hiparnt.ihpr2[11]=0;
  //--------------------------------------------------------------
  hijset_(&EFRM, FRAME[FRAMEID], PROJ[PROJID], TARG[TARGID], &IAP, &IZP, &IAT, &IZT);
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  if(MODE == 0){
    std::string seedStr = std::to_string(randomNumber);
    std::string outputFileName = "output_" + seedStr + ".hepmc";
    std::ofstream outputStream(outputFileName);
    WriterAscii outputWriter(outputStream);
    //------------------------------------------------------------
    //------------------------------------------------------------
    std::stringstream IonP;
    IonP<< std::setw(3) << std::setfill('0') << IZP;
    IonP<< std::setw(3) << std::setfill('0') << IAP;
    std::string IonPP = IonP.str();
    IonPP = "100" + IonPP + "0";
    int IonBeamP = std::stoi(IonPP);
    //--------------------------------------------------------------
    std::stringstream IonT;
    IonT<< std::setw(3) << std::setfill('0') << IZT;
    IonT<< std::setw(3) << std::setfill('0') << IAT;
    std::string IonPT = IonT.str();
    IonPT = "100" + IonPT + "0";
    int IonBeamT = std::stoi(IonPT);
    //--------------------------------------------------------------
    int Nev_tot = 0;
    for (int IE = 0; IE <= N_EVENT; ++IE) {    
      //------------------------------------------------------------
      hijing_(FRAME[FRAMEID], &BMIN, &BMAX);
      //------------------------------------------------------------
      if(himain1.natt <= 1) continue;
      if(hiparnt.hint1[18] > 1.2*(RIonP + RIonT) ) continue;
      //------------------------------------------------------------
      GenEvent Hep_evtent;
      Hep_evtent.set_units(Units::GEV,Units::MM);
      Hep_evtent.set_event_number(Nev_tot);
      //------------------------------------------------------------
      //1 fm = 1.0e-12 mm
      float ImbactP = hiparnt.hint1[18]*FMtoMM;
      //------------------------------------------------------------
      HepMC3::GenParticlePtr  BeamPP= std::make_shared<GenParticle>( FourVector( 0.0, 0.0,  EFRM/2., EFRM/2.), IonBeamP, 4);
      BeamPP->set_generated_mass(1.834732e+02);
      
      HepMC3::GenParticlePtr  BeamPT= std::make_shared<GenParticle>( FourVector( 0.0, 0.0, -EFRM/2., EFRM/2.), IonBeamT, 4);
      BeamPT->set_generated_mass(1.834732e+02);
      
      Hep_evtent.set_beam_particles(BeamPP, BeamPT);
      //------------------------------------------------------------
      int   npart0     = 0;
      float vxpart[600]={0};
      float vypart[600]={0};
      
      HepMC3::GenVertexPtr BeamP_vertice = std::make_shared<GenVertex>(FourVector( 0.5*ImbactP,0.0,0.0,0.0));
      BeamP_vertice->add_particle_in(BeamPP);
      for(int ip=0; ip<IAP; ip++){
	int idp = 0;
	if(ip < IZP) idp = 2212;
	else         idp = 2112;
	BeamP_vertice->add_particle_out( std::make_shared<GenParticle>( FourVector( histrng.pp[0][ip],
										    histrng.pp[1][ip],
										    histrng.pp[2][ip],
										    histrng.pp[3][ip]),
									idp, 0 ) );
	if(histrng.nfp[4][ip] > 1){
	  vxpart[npart0] = hijcrdn.yp[ip][0];
	  vypart[npart0] = hijcrdn.yp[ip][1];
	  npart0++;
	}
	
      }
      Hep_evtent.add_vertex(BeamP_vertice);
      //-----------------------------------------------------------
      HepMC3::GenVertexPtr BeamT_vertice = std::make_shared<GenVertex>(FourVector(-0.5*ImbactP,0.0,0.0,0.0));
      BeamT_vertice->add_particle_in (BeamPT);
      for(int it=0; it<IAT; it++){
	int idt = 0;
	if(it < IZT) idt = 2212;
	else         idt = 2112;
	BeamT_vertice->add_particle_out( std::make_shared<GenParticle>( FourVector( histrng.pt[0][it],
										    histrng.pt[1][it],
										    histrng.pt[2][it],
										    histrng.pt[3][it]),
									idt, 0 ) );
	if(histrng.nft[4][it] > 1){
	  vxpart[npart0] = hijcrdn.yt[it][0];
	  vypart[npart0] = hijcrdn.yt[it][1];
	  npart0++;
	}      
	
      }
      Hep_evtent.add_vertex(BeamT_vertice);
      //------------------------------------------------------------
      //cout<<hijcrdn.yt[it][0]<<" "<<hijcrdn.yt[it][1]<<" "<<hijcrdn.yt[it][2]<<"\n";
      //cout<<histrng.pp[0][it]<<" "<<histrng.pp[1][it]<<" "<<histrng.pp[2][it]<<" "<<histrng.pp[3][it]<<"\n";
      //------------------------------------------------------------
      int   Npart = 0;
      float Mx = 0.0,My=0.0;
      for(int i=0; i< (int) npart0; i++){
	Mx    += vxpart[i];
	My    += vypart[i];
	Npart ++;
      }
      if(Npart == 0) continue;
      Mx /=((float) Npart);
      My /=((float) Npart);
      
      float Rn[10],Xn[10],Yn[10];
      for(int i=0; i< 10; i++){
	Rn[i]=0.; Xn[i]=0.;Yn[i]=0.;
      }
      
      for(int i=0; i< (int) npart0; i++){
	float Vxd   = vxpart[i] - Mx;
	float Vyd   = vypart[i] - My;
	
	float r     = sqrt( pow(Vxd,2.) + pow(Vyd,2.) );
	float p     = atan2(Vyd , Vxd);
	
	for(int ni=0; ni<10; ni++){
	  float h = (float)1.0*ni + 1.0;
	  float w;
	  if(ni==0) w = 3.0;
	  else      w = h;
	  Rn[ni] += pow(r,w);
	  Xn[ni] += pow(r,w)*cos(h*p);
	  Yn[ni] += pow(r,w)*sin(h*p);
	}
      }
      
      double Psi[10]={0};
      double Ecc[10]={0};
      for(int ni=0; ni<10; ni++){
	float h  = (float)1.0*ni + 1.0;
	Xn[ni]  /=((float) Npart);
	Yn[ni]  /=((float) Npart);
	Rn[ni]  /=((float) Npart);
	
	Psi[ni] = (atan2(Yn[ni],Xn[ni]) + Pi )/h;
	Ecc[ni] = sqrt(pow(Xn[ni],2.) + pow(Yn[ni],2.))/Rn[ni];
      }
      //------------------------------------------------------------
      //------------------------------------------------------------
      //------------------------------------------------------------
      //------------------------------------------------------------
      vector<HepMC3::GenParticlePtr> Hep_particles;
      vector<int>                    Hep_mother_ids;
      vector<HepMC3::GenVertexPtr>   Hep_prods;
      
      Hep_particles.clear();
      Hep_mother_ids.clear();
      Hep_prods.clear();
      //------------------------------------------------------------
      //------------------------------------------------------------
      for(int part=0; part<himain1.natt; part++){
	//---------------------------------------------------------
	int mid   = himain2.katt[2][part] - 1;
	//---------------------------------------------------------
	int    id = himain2.katt[0][part];
	
	double px = himain2.patt[0][part];
	double py = himain2.patt[1][part];
	double pz = himain2.patt[2][part];
	double en = himain2.patt[3][part];
	
	double vx = himain3.vatt[0][part]*FMtoMM;
	double vy = himain3.vatt[1][part]*FMtoMM;
	double vz = himain3.vatt[2][part]*FMtoMM;
	double tt = himain3.vatt[3][part]*FMtoMM;
	//---------------------------------------------------------
	int status = himain2.katt[3][part];
	if(status <=10 && status>0 ) status = 1;
	if(status <=20 && status>10) status = 2;
	//---------------------------------------------------------
	Hep_particles.push_back(std::make_shared<GenParticle>( FourVector( px, py, pz, en), id, status));
	Hep_prods.push_back(std::make_shared<GenVertex>(FourVector(vx,vy,vz,tt)));
	Hep_mother_ids.push_back(mid);
      }
      //------------------------------------------------------------
      //------------------------------------------------------------
      for (unsigned int ipart = 0; ipart<Hep_particles.size(); ipart++) {
	//---------------------------------------------------------
	int mid = Hep_mother_ids[ipart];
	//---------------------------------------------------------
	double pT = sqrt(Hep_particles[ipart]->momentum().px()*Hep_particles[ipart]->momentum().px() +
			 Hep_particles[ipart]->momentum().py()*Hep_particles[ipart]->momentum().py());
	if(Hep_particles[ipart]->status()==1 && pT==0.0) continue;
	//---------------------------------------------------------
	if(mid < 0 ){
	  Hep_evtent.add_particle(Hep_particles[ipart]);
	}else{
	  GenVertexPtr Prod_vertice = std::make_shared<GenVertex>();
	  Prod_vertice = Hep_prods[mid];
	  Prod_vertice->add_particle_in (Hep_particles[mid]  );
	  Prod_vertice->add_particle_out(Hep_particles[ipart]);
	  Hep_evtent.add_vertex(Prod_vertice);
	}
	//--------------------------------------------------------
      }
      //------------------------------------------------------------
      //------------------------------------------------------------
      double Sigma_NN = hiparnt.hint1[11]*pow(10,9);//In the unit of Pb
      //------------------------------------------------------------
      Hep_evtent.weights().push_back(1.0);
      //------------------------------------------------------------
      float cent = 0.0;
      cent       = 100*pow((hiparnt.hint1[18]/(2.*RPT)),2.0);
      //------------------------------------------------------------
      std::shared_ptr<GenHeavyIon> heavy_ion = std::make_shared<GenHeavyIon>();
      Hep_evtent.add_attribute("GenHeavyIon",heavy_ion);
      heavy_ion->set(
		     himain1.jatt,
		     himain1.np,
		     himain1.nt,
		     himain1.n0+himain1.n01+himain1.n10+himain1.n11,
		     0,
		     0,
		     himain1.n01,
		     himain1.n10,
		     himain1.n11,
		     hiparnt.hint1[18],
		     0.0,
		     0.0,
		     Sigma_NN,
		     cent    );
    
      for(int i=0; i<10; i++){
	heavy_ion->participant_plane_angles[i] = Psi[i];
	heavy_ion->eccentricities[i]           = Ecc[i];
      }
      Hep_evtent.set_heavy_ion(heavy_ion);
      //---------------------------------------------------------
      std::shared_ptr<GenCrossSection> cross_section = std::make_shared<GenCrossSection>();
      Hep_evtent.add_attribute("GenCrossSection",cross_section);
      cross_section->set_cross_section(Sigma_NN, 0.000000001*Sigma_NN);
      //---------------------------------------------------------
      outputWriter.write_event(Hep_evtent);
      Hep_evtent.clear();
      //---------------------------------------------------------
      //---------------------------------------------------------
      Nev_tot ++ ;
      //---------------------------------------------------------
    }//end of event loop
    //-----------------------------------------------------------
    outputWriter.close();
    //-----------------------------------------------------------
  }
  //-------------------------------------------------------------
  else{
    //-----------------------------------------------------------
    int    Nnpar     = 0;
    int    Nmult     = 0;
    float  Impac     = 0.;
    float  Psin[10]  ={0};
    float  Eccn[10]  ={0};

    int    fId[80000]= {0};
    float  fPx[80000]= {0};
    float  fPy[80000]= {0};
    float  fPz[80000]= {0};
    float  fEn[80000]= {0};
    
    char ROOTName[500];
    sprintf(ROOTName,"output_%d.root",randomNumber);
    TFile f(ROOTName,"recreate");
        
    TTree *tree = new TTree("tree","tree");
    
    tree->Branch("Nnpar" ,&Nnpar,"Nnpar/I");
    tree->Branch("Nmult" ,&Nmult,"Nmult/I");
    tree->Branch("Impac" ,&Impac,"Impac/F");
    //------------------------------------------------------------
    tree->Branch("Psin",Psin,"Psin[10]/F");
    tree->Branch("Eccn",Eccn,"Eccn[10]/F");
    //------------------------------------------------------------
    tree->Branch("fId",fId,"fId[Nmult]/I");
    tree->Branch("fPx",fPx,"fPx[Nmult]/F");
    tree->Branch("fPy",fPy,"fPy[Nmult]/F");
    tree->Branch("fPz",fPz,"fPz[Nmult]/F");
    tree->Branch("fEn",fEn,"fEn[Nmult]/F");
    //------------------------------------------------------------
    //------------------------------------------------------------
    int Nev_tot = 0;
    for (int IE = 0; IE <= N_EVENT; ++IE) {    
      //----------------------------------------------------------
      hijing_(FRAME[FRAMEID], &BMIN, &BMAX);
      //----------------------------------------------------------
      if(himain1.natt <= 1) continue;
      if(hiparnt.hint1[18] > 1.2*(RIonP + RIonT) ) continue;
      //----------------------------------------------------------
      int   npart0     = 0;
      float vxpart[600]={0};
      float vypart[600]={0};
      
      for(int ip=0; ip<IAP; ip++){
	if(histrng.nfp[4][ip] > 1){
	  vxpart[npart0] = hijcrdn.yp[ip][0];
	  vypart[npart0] = hijcrdn.yp[ip][1];
	  npart0++;
	}
      }
      for(int it=0; it<IAT; it++){
	if(histrng.nft[4][it] > 1){
	  vxpart[npart0] = hijcrdn.yt[it][0];
	  vypart[npart0] = hijcrdn.yt[it][1];
	  npart0++;
	}      
      }
      //------------------------------------------------------------
      int   Npart = 0;
      float Mx = 0.0,My=0.0;
      for(int i=0; i< (int) npart0; i++){
	Mx    += vxpart[i];
	My    += vypart[i];
	Npart ++;
      }
      if(Npart == 0) continue;
      Mx /=((float) Npart);
      My /=((float) Npart);
    
      float Rn[10],Xn[10],Yn[10];
      for(int i=0; i< 10; i++){
	Rn[i]=0.; Xn[i]=0.;Yn[i]=0.;
      }
    
      for(int i=0; i< (int) npart0; i++){
	float Vxd   = vxpart[i] - Mx;
	float Vyd   = vypart[i] - My;
      
	float r     = sqrt( pow(Vxd,2.) + pow(Vyd,2.) );
	float p     = atan2(Vyd , Vxd);
      
	for(int ni=0; ni<10; ni++){
	  float h = (float)1.0*ni + 1.0;
	  float w;
	  if(ni==0) w = 3.0;
	  else      w = h;
	  Rn[ni] += pow(r,w);
	  Xn[ni] += pow(r,w)*cos(h*p);
	  Yn[ni] += pow(r,w)*sin(h*p);
	}
      }
    
      double Psi[10]={0};
      double Ecc[10]={0};
      for(int ni=0; ni<10; ni++){
	float h  = (float)1.0*ni + 1.0;
	Xn[ni]  /=((float) Npart);
	Yn[ni]  /=((float) Npart);
	Rn[ni]  /=((float) Npart);
      
	Psi[ni] = (atan2(Yn[ni],Xn[ni]) + Pi )/h;
	Ecc[ni] = sqrt(pow(Xn[ni],2.) + pow(Yn[ni],2.))/Rn[ni];
      }
      //----------------------------------------------------------
      //----------------------------------------------------------
      int Nch = 0;
      //------------------------------------------------------------
      for(int part=0; part<himain1.natt; part++){
	//---------------------------------------------------------
	int    id = himain2.katt[0][part];
	double px = himain2.patt[0][part];
	double py = himain2.patt[1][part];
	double pz = himain2.patt[2][part];
	double en = himain2.patt[3][part];
	//---------------------------------------------------------
	int status = himain2.katt[3][part];
	if(status <=10 && status>0 ) status = 1;
	if(status <=20 && status>10) status = 2;
	//---------------------------------------------------------
	if(status != 1) continue;
	double pt = sqrt(pow(px,2) + pow(py,2) );
	if(status == 1 && pt == 0.0) continue;
	//---------------------------------------------------------
	fId[Nch] = id;
	fPx[Nch] = px;
	fPy[Nch] = py;
	fPz[Nch] = pz;
	fEn[Nch] = en;
	Nch ++;
      }
      //----------------------------------------------------------
      if(Nch > 1){
	Nmult = Nch;
	Nnpar = npart0;
	Impac = hiparnt.hint1[18];
	
	for(int ni=0; ni<10; ni++){
	  Psin[ni] = Psi[ni];
	  Eccn[ni] = Ecc[ni];
	}
	tree->Fill();
      }

      Nev_tot ++ ;
    }
    //------------------------------------------------------------
    tree->Write();
    f.Close();
  }
  //-------------------------------------------------------------
  //-------------------------------------------------------------
  std::cout<<"The end of the code"<<"\n";
  //-------------------------------------------------------------
  //-------------------------------------------------------------  
  return 0;
}
