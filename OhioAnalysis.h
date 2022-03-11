#include "kinfit.h"
#include "sstream"
#include <algorithm>
#include "TH1F.h"
#include "TH2F.h"
#include "TBrowser.h"
#include "TMinuit.h"
#include "TEnv.h"
#include "Math/Minimizer.h"
#include "Math/GSLMinimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include <iostream>
#include <string>
#include "TCutG.h"
#include "TList.h"
#include "TFile.h"
#include "TCanvas.h"
#include <vector>
#include "TF1.h"
#include "TSpectrum.h"
#include "TLegend.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TH3F.h"
#include "fstream"
#include "Particle.h"
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */
double pi=4.*atan(1.);
using namespace std;
class Analysis {
public:
	Analysis(bool);
	~Analysis();
	void DefineHists();
	void Execute(TString);
	void WriteHist(TString);
	void DrawTheory();
	void ReadMap();
	void Dist(int,int,int,double&,double&,double&);
	void DistNoSmear(int,int,int,double&,double&,double&);
	void SetResponse();
	void SetResponse2();
	void SetResponse3();
	void SetExperiment(int);
	void HoughTransform(int);
	void FindPlane();
	void FillRunRange();
	void FitDecay();
	void SetRunNumber(int);
	int runnumber;
	int HoughTransformEvent();
	double GetPressure(double);
	void HoughTransform3D();
	double FitLine();
	int AnalyzeEvent(int);
	void HoughTransform2();
	const int nsamples=10;
	void HoughTransformGeneric(TH2F*,double&,double&);
//	double FitSaturatedWaveform(double,TH1F*,int);
	void StripsChains();
	void CurvedArm(int);
	void DefineEventHistograms();
	void GetCenter(int,int,int,int&,int&,int&,double);
	void GetEnd(double&,double&,double&);
	void GetArm(int,int,int,int&,int&,int&);
	int CleanY(int&,int&,int&,int&,int&,int&,int,int,int,bool&);
	void RemoveOutliers(int,int,int);
	void SetHistName(TString);
	double GetRange(double);
	double GetRange12C(double);
	double GetRange16O(double);
	double GetEnergy(double);
	double GetEnergy12N(double);
	double GetEnergy4He(double);
	double GetEnergy9Be(double);
	double GetEnergy13C(double);
	double GetEnergy12C(double);
	double * Get23Arm(double arr[3][2],int,int,int,int,int,int,bool);
	void GetPlane();
	void Ransac();
	void Ransac2();
	void Ransac3();
//	void CalcChisqrdDist(Int_t &, Double_t *, Double_t &, Double_t *, Int_t);
	double CalcChisqrdDist(const double *);
	double CalcChisqrdCurve(const double *);
	double CalcChisqrdLine(const double *);
	TFile *filetrack;
	double kinfitchisqrd;
	TString GetHistName();
	void FindSecondThirdArm();
	void FindBeamStop();
	void OtherArms();
	void GetWiggliness();
	void FillEnergyRange();
	void ReadCalibrationFile(TString);
	void AnalyzeSi();
	double RemoveBG();
	void EventIdentificationFunction();
	void Cross(double,double,double,double,double,double,double&,double&,double&);
	void EnergyTimeBeam();
	double GetWeight(double,double);
	double arm_point[3];
	double arm_pix1[3];
	double arm_pix2[3];
	double arm_point1[3];
	double arm_point2[3];
	double arm_point3[3];
	double scattering_theta[3];
	double arm_energy[3];
	double chi_linear,chi_triple;
	double endpoint[2][3];
	double GetLabAngle12C(double);
	double GetLabAngle16O(double);
	double GetLabAngle16Op(double);
	double GetCMAngle12C(double);
	double GetCMAngle16O(double);
	double GetCMAngle16Op(double);
	double p_vertex[3]={0};
	double p_end1[3]={0};
	double p_end2[3]={0};
	double p_chisqrd=0.;
	void KinematicFitting();
	void myfunc();
	void FitSaturated(double&,double&,double,TH1F*);
	Double_t pulseshape(Double_t[],Double_t[]);
	void HalfEvent(int);
	void Deconvolution();
	void GetDot(double,double,double,double,double,double);
	void GetOverlapPoint(double&,double&,double&);
	void ColinearDecayTest();
	void Fits();
	void Flatten(int,int,int);
//	double CalcChisqrdDist(const double *);
	void FancyCenter();
	void NearestNeighbour();
	void DirectDalitz();
	void FindReactionPlane2D();
	double vertexy=0.;
	double Fitnalpha();
	double Eleft=0,Eright=0.;
	bool verbose=false;
	double energy[3];
	double chisqrd_nalpha=0;
	double chisqrd_nalpha_theta=0;
	double chisqrd_nalpha_dist=0;
	int select_nalpha=1;
	int looper=0;
	int globaleventno=0;
	int ExecuteList(TString);
	int centralhits[6][128],centralhitsT[6][128];
	int centralhitsevt[6][128],centralhitsTevt[6][128];
	int starting_event,draw_events,drawing,max_draw,timetemp,good_evt,doublecount;
	int map_x[4][4][64],map_y[4][4][64];
	int simap_x[1][1][64],simap_y[1][1][64],EventNo;
	int evtno,waveformno,deconvno,deconvnob,event_selection;
	int maxx,maxy,maxz,step_count;
	static const int c_events=100;
	double /*driftv=0.0075,*/totalEnergy=0,allchanE=0,allchanE2=0.;
	double pressure=50.;//0.0075 mm/ns for 900 V across 13 cm 50 Torr CO2
//	double pressure=100.;//0.0075 mm/ns for 900 V across 13 cm 50 Torr CO2
	double driftv=0.0075;//fudge factor for exp 50 Torr
//	double driftv=0.0043;//300 bin spread seen in data
	double mc_thresh,E_thresh,E_threshb,dist_diff;
	double L1Agradient,Ex;
	bool good,fulleventreadout,mixedevent;
	float xgoodpeaks[100];
	const int striparsize=140;//make this array large enough to fit the width of the pseudopixels (140)
	double mmstripsE[140][2][2],mmchainsE[140][2],mmstripsT[140][2][2],mmchainsT[140][2],mmstripsE_2[140][2][2],mmchainsE_2[140][2],mmstripsW[140][2][2],mmchainsW[140][2];
	double length[3]={0.};
	double arm1_E=0,arm2_E=0,arm3_E=0;
	static const int maxsize=3264;
	bool debug,takeall,fpnadjust,halfevent;
	double chisqrdHoyle=0;
	double linetheta=0.;
	double dE_light,dE_heavy;
	double dE_light_[3],dE_heavy_[3];
	double linedist=0.;
        double ftheta[2][2]={{0.,0.},{0.,0.}};//for two light/heavy combos for theta/psi

        Int_t mmWaveformY[maxsize][512];//array to hold waveform
        Int_t CorWaveform[maxsize][512];//array to hold waveform
        Int_t mmPixX[160];//pixel hit
        Int_t mmPixY[160];//pixel hit
        Int_t mmCobo[maxsize],mmAsad[maxsize],mmAget[maxsize],mmChan[maxsize],mmDecayNo[maxsize],mmHit,mmFrameNo[maxsize],maxval[maxsize],mmEventIdx;
        Int_t mmMul,evt,events;
        Int_t mmSCA1,mmSCA2,mmSCA3,mmSCA4,mmSCA5,mmSCA6;
	Double_t mmD2PTime;
	double FullRange[200];
	double eventchisqrd,dheight,maxheight,minheight;
	double mmenergy[1024][2],mmtime[1024][2],mmenergy_2[1024][2],mmwidth[1024][2];
	double FPNwaveform[4][4][256]={{{0.}}};
	double padenergy,peakEnergy=0,peaking,peaktime=0,sumE;
	bool isFPN[maxsize],domore=false;
	double fit_theta[3][2];
	double theory_theta[2];
	double wave[512];
	double beamwave[256];
        double responseshift[512];
	double rlight[3];
	double rheavy[3];
	double rtlight[3];
	double rtheavy[3];
        int offset2=0;
        int source_offset=0;
	int good_double=0;
	int beamx=0,beamy=0,beamz=0;
	double dmaxx,dmaxy,dmaxz,scatterangle;
	double scaletime=100*driftv;//100ns per bin * driftV = distance
	double maxdist=0,Ebeam=0.,firstpadE[6],secondpadE[6],maxdist2=0.,maxdist3=0.;
	double baseline_set=50.;
	double totalEnergy_glob=0,ICenergy=0.,Elength;
	double E_scaling,decay_x,decay_y,decay_z;
	double decay_x2,decay_y2,decay_z2;
	double theta2decay,theta3decay,spread,highy,lowy;
	double arm1_x,arm1_y,arm1_z;
	double arm2_x,arm2_y,arm2_z;
	double arm3_x,arm3_y,arm3_z;
	double darm1_x,darm1_y,darm1_z;
	double darm2_x,darm2_y,darm2_z;
	double darm3_x,darm3_y,darm3_z;
	double palpha[3][3],palphakf[3][3],excitationenergy;
	double Ealpha[3],Ralpha[3],Salpha[3],Ralphakf[3],Ealphakf[3];
	double thetaalpha[3],phialpha[3],errthetaalpha[3],errphialpha[3],Ptalpha[3],errPtalpha[3];
	double centr_dist_x,centr_dist_y,centr_dist_z,beam_dist_x,beam_dist_y,beam_dist_z,centr_dist_x2,centr_dist_y2,centr_dist_z2;
	int fexperiment=0,largest_time=0;
	int remove[3]={0};
	double HoughLine[3][3],Normal[3],neutronE;
	double gain[20],offset[20];
	bool canstripmatchL,canstripmatchR,saturated;
        int mult_decay=0,saturatedchnls=0;
        int firstE=0,secondE=0,second_E=0;
        int largesty=0,largesty_b=0,oddframes=0,evenframes=0;
	double beam_divergence,decay_divergence,totalbeamE=0.,HoughCenterX,HoughCenterY,HoughCenterZ,integratedE;
	double dtheta[3];
	double dlength[3][2];
	double integratewf;
	TString outputdefaultname;
	TTree *tree;
	TTree *treetrack;
	TH1I *hevents;
	TH1I *hcount;
	TH1F *htheta12C,*htheta16O,*htheta12CCM,*htheta16OCM,*htheta16Op,*htheta16OpCM;
	TH1F *h1,*h2,*h2_g,*h3,*h3_g,*h3_g0,*h3_g1,*h3_g2,*h4,*h5,*hE_all,*hpulse_chn[c_events]={NULL},*hpulse_chndc[c_events]={NULL},*hdedx[c_events]={NULL},*hbg,*htime,*htime2,*hmaxtrack,*hmaxtrack_g,*hmaxtrackfit,*hmaxtrackkf,*hmaxtracktheory,*hdoublecandidate[1000],*hpeaks,*hmult,*hrndwalk,*hrndwalkx,*hrndwalky,*hdot[c_events],*hangle=0,*hlength_all=0,*htotal_length=0,*htotal_lengthkf,*htotal_lengththeory,*hEx,*hExbest,*hcurdedx=0,*hextrapolate,*hmaxarc_length;
	TH2F *htrackxysmooth[c_events]={NULL},*htrackxy[c_events]={NULL},*htrackyz[c_events]={NULL},*htrackxz[c_events]={NULL},*hpulse[c_events]={NULL},*hpulse_b[c_events]={NULL},*hHits[c_events]={NULL},*hHits_side[c_events]={NULL},*hEtime,*hallhits,*hallhits2,*hlengthvsE,*hwaveform[c_events]={NULL},*hwaveform_dc[c_events]={NULL},*hwaveform_hit[c_events],*hwaveformb[c_events],*hwaveform_dcb[c_events],*hwaveform_hitb[c_events],*htrackxyvert[c_events]={NULL},*htrackyzvert[c_events];
	TH2F *hbeamtime,*hbeam[c_events],*hbeam2[c_events],*hbeamstop,*htheta2theta3,*htheta2theta3fit,*htheta2theta3kf,*htheta2theta3theory,*hDalitz,*hDalitztheory,*hDalitzkf,*hthetapsi,*hbeamEdE,*hEtotntot,*hmaxtrack_height,*harm_E2D,*harm_height,*hHoughXY2[c_events]={NULL},*hHoughYZ2,*hHoughXZ2[c_events]={NULL},*hHoughXY[c_events]={NULL},*hHoughYZ[c_events]={NULL},*hHoughXZ[c_events]={NULL};
	TH2F *htrackcor,*hangle_arms[c_events],*hpeaktotal,*hEx1Ex2,*hEx_d,*hEx_evt,*hsuper[c_events],*hdEE,*hdeL,*hcorrelation,*hdeltat[c_events]={NULL},*hflat[c_events]={NULL},*harm1_xy[c_events]={NULL},*harm1_yz[c_events],*hfitbeamXY[c_events]={NULL},*hfitbeamYZ[c_events]={NULL};
	TH2F *hdistxy[c_events]={NULL},*hdistyz[c_events]={NULL},*hdistxz[c_events]={NULL},*hdistxy1[c_events]={NULL},*hdistxy2[c_events]={NULL},*hdistxy2a[c_events]={NULL},*hdistyz2[c_events]={NULL},*hdistxz2[c_events]={NULL},*hdistxy2b[c_events],*hdistxy2c[c_events],*hdistyz2c[c_events],*hdistxz2c[c_events],*hdistxy3[c_events]={NULL},*hdistyz1[c_events]={NULL},*hdistyz2a[c_events]={NULL},*hdistyz2b[c_events],*hdistyz3[c_events]={NULL},*hdistxz1[c_events]={NULL},*hdistxz2a[c_events]={NULL},*hdistxz2b[c_events],*hdistxz3[c_events]={NULL},*hstage1[c_events]={NULL},*hstage2[c_events]={NULL},*hstage3[c_events]={NULL},*hbeam_c,*hflat_all;
	TH3F *htrack3D[c_events]={NULL},*hcurrenttrack3D=0,*hcurrenttrack3D_2=0,*hcurrenttrack3D_3=0,*hHoughXYZ=0,*hHoughPlane1,*hHoughPlane2,*hbeamtrack,*htrackarm1,*htrackarm2,*htrackarm3;
	TH2F *harmxy[c_events],*harmyz[c_events],*harmxz[c_events];
	TH1F *henergy[c_events]={NULL},*hbeamenergy,*harm_E,*hD2Ptime,*hnearestneighbour[c_events]={NULL},*hmults[c_events]={NULL},*hHoughprojXY[c_events]={NULL},*hHoughprojYZ[c_events]={NULL},*hHoughprojXZ[c_events]={NULL},*hchainscheckL,*hchainscheckR,*hchainsevt[c_events]={NULL},*hwave,*hwave2,*hIC;
	TH2F *hmaxchains[c_events]={NULL},*hmaxstrips[c_events]={NULL},*hbeamEl,*hHoughXYb,*hHoughXZb,*hHoughYZb;
	TH1F *hbdivxy,*hbdivyz,*hBe,*hBekf,*hEtotal,*hpx,*hpy,*hpz,*hpx_kf,*hpy_kf,*hpz_kf,*hEventID,*hEventIDcut,*htheta2,*htheta2kf,*hfitchisqrd;
	TH2F *hdriftcheck,*hdist12,*hExl,*hrndwalktime,*hrndwalk2D;
	TH2F *hdiv,*hallsihits,*hDSSD,*hEventID2D,*hEventIDheight,*hEventIDposx,*hEventIDlength,*h5l,*h5t,*h5h,*hmultspread,*hmultspread2,*hmultE,*hmultE2;
	TH1F *hEventIDgoodevent,*hcolinear,*hbeamstopx,*hbeamstopy,*hdisplacement,*hneutronenergy;
	TH2F *hHoylestopY,*hd2ptimebeamy,*hE_ally,*hEtheta,*hEthetatheory,*hneutronenergy2D;
	TH1F *hSilicon[20],*hSilicong[20],*hSilicongg[20],*hChi1,*hChi2,*hwidth,*hQvalue12C,*hQvalue16O;
	TH2F *hrangeE1,*hrangeE2,*hrangeE3,*hrange12,*hrange23,*hkfxy[c_events]={NULL},*hkfyz[c_events]={NULL},*hkfxz[c_events]={NULL},*hWiggle1,*hWiggle2,*hWiggle3,*hWiggleL1,*hWiggleL2,*hWiggleL3,*hWiggleevt[c_events]={NULL},*hcurve,*hRtotRmax,*hscatterEl;
	TH2F *hErEs,*hE1E2,*hcurrentYZ,*hlinearity,*hthetatheta,*hthetathetaraw,*hthetatheta12C,*hthetatheta16O,*hthetatheta16Op,*hthetathetatheory,*hQvalue2D,*hEavg,*hchisqrd,*hthetathetatheoryCM;
	TH1F *hchisqrdth12C,*hchisqrdth16O,*hchisqrdth16Op;
	TH1F *hchisqrdE12C,*hchisqrdE16O,*hchisqrdE16Op;
	TH1F *hselect,*hselect12C16O,*hselect12C16Op,*hselect16O16Op,*hconf12C,*hconf16O,*hconf16Op,*hchisqrdnalpha,*hconf212C,*hconf216O,*hconf216Op;
	TH2F *hselect2,*hselect3,*hconf12C2D,*hconf16O2D,*hconf16Op2D;
	TH2F *hchi12Ctheta,*hchi16Otheta,*hchi16Optheta,*helscat;
	TH2F *hlightR,*hheavyR;
        bool readmc,twopmode;
	double response[512];
        TRandom3 *rnd3;
	TF1 *fit,*mygaussian;
	TCutG* gcut_1;
	bool colinearflag=false;
	double colinearityconfidence;
	static const int ERentries=46;
	double EnergyRange[2][ERentries];
	int ingate=0;
	double d2r=atan(1.)/45.;//4*atan(1.)=pi-->pi/180=atan(1.)/45
	vector<double> pointx;
	vector<double> pointy;
	vector<double> pointz;
	vector<double> pointE;
	vector<double> pointx2;
	vector<double> pointy2;
	vector<double> pointz2;
	vector<double> pointE2;
	vector<double> theta12C;
	vector<double> psi12C;
	vector<double> theta16O;
	vector<double> psi16O;
	vector<double> E16O_light;
	vector<double> E16O_heavy;
	vector<double> theta16Op;
	vector<double> psi16Op;
	vector<double> E16Op_light;
	vector<double> E16Op_heavy;
	vector<double> E12C_light;
	vector<double> E12C_heavy;
};

Analysis::Analysis(bool readmc_) {
/////////////////Set up events to read/events to draw e.t.c.//////////////////////////////////////////////////////////////////////
//	events=120;
//	mygaussian = new TF1("mygaussian","[0]*TMath::Exp(-TMath::Power(x-[1],2)/(2.*[2]*[2]))",0,512);
	mygaussian = new TF1("mygaussian","[0]*TMath::Exp(-TMath::Power(x-[1],2)/(2.*[2]*[2]))+[3]+x*[4]",0,512);
	readmc=readmc_;
	rnd3 = new TRandom3();
	cout<<"New analysis object defined\n";
	max_draw = c_events;
	starting_event=0;
	draw_events = max_draw;
	drawing=0;
	timetemp=256;
	if(readmc)timetemp=512;
//	scaletime=1.1;
//	if(readmc)scaletime=0.27;
	if(readmc)scaletime=0.27*2.;
	good_evt=0;
	hevents = new TH1I("hevents","Event reduction",20,0.,20);
	hcount = new TH1I("hcount","Count - fExperiment",20,0.,20);
	h1 = new TH1F("h1","Total energy",1000,0.,1e7);
	hpeaktotal = new TH2F("hpeaktotal","Peak energy against total energy",100,0.,3e5,100,0.,3e6);
	h2 = new TH1F("h2","Peak energy (MeV)",1000,0.,1e6);
	h2_g = new TH1F("h2_g","Peak energy (MeV) gated",150,5.,20.);
	h3 = new TH1F("h3","Total energy",1000,0.,2e6);
	h3_g = new TH1F("h3_g","Total energy (MeV) gated on track length",1000,0.,2e6);
	h3_g0 = new TH1F("h3_g0","Total energy (MeV) gated on 12C",1000,0.,2e6);
	h3_g1 = new TH1F("h3_g1","Total energy (MeV) gated on 16O",1000,0.,2e6);
	h3_g2 = new TH1F("h3_g2","Total energy (MeV) gated on 16Op",1000,0.,2e6);
	h4 = new TH1F("h4","Total integrated energy all channels",200,0,2e6);
	h5 = new TH1F("h5","Weighted average of integrated and peak energy spectrum",200,6.,10.);
	h5l = new TH2F("h5l","Weighted average of integrated and peak energy spectrum vs max track length",200,6.,8.,200,0.,50);
	h5t = new TH2F("h5t","Weighted average of integrated and peak energy spectrum vs Ex from total track length",200,6.,8.,200,6.,8);
	h5h = new TH2F("h5h","Weighted average of integrated and peak energy spectrum vs height of track",200,6.,8.,150,0.,150);
	hE_all = new TH1F("hE_all","Peak energy (no biased gating)",150,5,20);
	hE_ally = new TH2F("hE_ally","Peak energy (no biased gating) vs beam y",150,5,20,128,0,128);
	hIC = new TH1F("hIC","Ion counter peak signal size",500,0,6000);
	hdEE = new TH2F("hdEE","Energy deposited in MM chains vs Si front",500,0.,10000,100.,0.,10.);
	hdeL = new TH2F("hdeL","Energy deposited in IC against the stopping point of beam",100,0.,8000,128,0.,400);
	hbeamenergy = new TH1F("hbeamenergy","Energy of the 12N beam from the range (window-stop)",100,0.,30.);
	hBe = new TH1F("hBe","8Be excitation energy",100,0,2.);
	hBekf = new TH1F("hBekf","8Be excitation energy kin. fit",100,0,2.);
	hEx_d = new TH2F("hEx_d","Excitation energy (MeV) against Y center",200,5.,10.,128,0,128);
	hEx_evt = new TH2F("hEx_evt","Excitation energy (MeV) against percentage through run",200,5.,15.,50,0,100);
	hEx1Ex2 = new TH2F("hEx1Ex2","Ex from total vs peak",200,5.,15.,200,5.,15.);
	hmult = new TH1F("hmul","Multiplicity plot",1024,0,1024);
	hmultspread = new TH2F("hmultspread","Multiplicity plot vs diff in Y",1024,0,1024,128,0,128);
	hmultspread2 = new TH2F("hmultspread2","Multiplicity plot vs diff in Y",1024,0,1024,128,0,128);
	hchisqrd = new TH2F("hchisqrd","12C (x) vs 16O (y)",500,0,100,500,0,100);
	hchisqrdth12C = new TH1F("hchisqrdth12C","12C chisqrd",1000,0,10);
	hchisqrdth16O = new TH1F("hchisqrdth16O","16O chisqrd",1000,0,10);
	hchisqrdth16Op = new TH1F("hchisqrdth16Op","16Op chisqrd",1000,0,10);
	hchisqrdE12C = new TH1F("hchisqrdE12C","12C chisqrd",1000,0,10);
	hchisqrdE16O = new TH1F("hchisqrdE16O","16O chisqrd",1000,0,10);
	hchisqrdE16Op = new TH1F("hchisqrdE16Op","16Op chisqrd",1000,0,10);
	hselect = new TH1F("hselect","chi_c - chi_o / (chi_c + chi_o)",100,-1,1);
	hselect12C16O = new TH1F("hselect12C16O","chi_c - chi_o / (chi_c + chi_o)",100,-1,1);
	hselect12C16Op = new TH1F("hselect12C16Op","chi_c - chi_o / (chi_c + chi_o)",100,-1,1);
	hselect16O16Op = new TH1F("hselect16O16Op","chi_c - chi_o / (chi_c + chi_o)",100,-1,1);
	hconf12C = new TH1F("hconf12C","chi_c - chi_o / (chi_c + chi_o)",100,-1,1);
	hconf16O = new TH1F("hconf16O","chi_c - chi_o / (chi_c + chi_o)",100,-1,1);
	hconf16Op = new TH1F("hconf16Op","chi_c - chi_o / (chi_c + chi_o)",100,-1,1);
	hconf212C = new TH1F("hconf212C","chi_c - chi_o / (chi_c + chi_o)",100,-1,1);
	hconf216O = new TH1F("hconf216O","chi_c - chi_o / (chi_c + chi_o)",100,-1,1);
	hconf216Op = new TH1F("hconf216Op","chi_c - chi_o / (chi_c + chi_o)",100,-1,1);
	hconf12C2D = new TH2F("hconf12C2D","chi_c - chi_o / (chi_c + chi_o)",100,-1,1,100,0,20);
	hconf16O2D = new TH2F("hconf16O2D","chi_c - chi_o / (chi_c + chi_o)",100,-1,1,100,0,20);
	hconf16Op2D = new TH2F("hconf16Op2D","chi_c - chi_o / (chi_c + chi_o)",100,-1,1,100,0,20);
	hselect2 = new TH2F("hselect2","chi_c - chi_o / (chi_c + chi_o) (x) vs min (y)",100,-1,1,100,0,50);
	hselect3 = new TH2F("hselect3","chi_c - chi_o / (chi_c + chi_o) (x) vs vertex y (y)",100,-1,1,128,0,224);
	hmultE = new TH2F("hmultE","Multiplicity plot vs Etot",1024,0,1024,1000,0,1e6);
	hmultE2 = new TH2F("hmultE2","Multiplicity plot vs Etot",1024,0,1024,1000,0,1e6);
	harm_E = new TH1F("harm_E","Total long arm energy",100,0.,.5);
	harm_E2D = new TH2F("harm_E2D","Total long arm energy vs length",100,0.,.5,100,0.,100.);
	htime = new TH1F("htime","D2P time",100,0,3e7);
	htime2 = new TH1F("htime2","D2P time (zoomed for beam+beam)",100,0,2e5);
	hbg = new TH1F("hbg","Background",timetemp,0,timetemp);
	hEtime = new TH2F("hEtime","Total energy vs d2p time",100,0.,5e4,100,0.,3e7);
	hallhits = new TH2F("hallhits","All hits in pixels",140,0,140,128,0,128);
	hallhits2 = new TH2F("hallhits2","All hits in position",240,0,245.,128,0,224.);
	hcorrelation = new TH2F("hcorrelation","Largest beam Y against largest decay Y",128,0,128,128,0,128);
	hmaxtrack = new TH1F("hmaxtrack","Max track length",256,0,256);
	hmaxtrack_g = new TH1F("hmaxtrack_g","Max track length gated on Hoyle",256,0,256);
	hmaxtrackfit = new TH1F("hmaxtrackfit","Max track length fitted",256,0,256);
	hmaxtrackkf = new TH1F("hmaxtrackkf","Max track length kinematic fitting",256,0,256);
	hmaxtracktheory = new TH1F("hmaxtracktheory","Max track length",256,0,256);
	hmaxtrack_height = new TH2F("hmaxtrack_height","Max track length versus the height traveled (to check drift vel)",256,0,256,100,0,100);
	hmaxarc_length = new TH1F("hmaxarc_length","Max track length using arc",256,0,256);
	hextrapolate = new TH1F("hextrapolate","Extrapolated arm length",256,0,256);
	htrackcor = new TH2F("htrackcor","Correction",128,0,128,128,0,128);
	hlength_all = new TH1F("hlength_all","All track lengths (mm)",256,0,256);
	htotal_length = new TH1F("htotal_length","Total track lengths (mm)",100,0,200);
	htotal_lengthkf = new TH1F("htotal_lengthkf","Total track lengths (mm)",100,0,200);
	htotal_lengththeory = new TH1F("htotal_lengththeory","Total track lengths (mm)",100,0,200);
	hlengthvsE = new TH2F("hlengthvsE","Max track length vs E",256,0,256,200,7.,14.);
	htheta12C = new TH1F("htheta12C","Angle for 12C",180,0,180);
	htheta16O = new TH1F("htheta16O","Angle for 16O",180,0,180);
	htheta16Op = new TH1F("htheta16Op","Angle for 16O(n,a1)",180,0,180);
	htheta12CCM = new TH1F("htheta12CCM","CM Angle for 12C",180,0,180);
	htheta16OCM = new TH1F("htheta16OCM","CM Angle for 16O",180,0,180);
	htheta16OpCM = new TH1F("htheta16OpCM","CM Angle for 16O",180,0,180);
	hEx = new TH1F("hEx","Excitation energy from arm lengths",1000.,0,10.);
	hExbest = new TH1F("hExbest","Excitation energy from arm lengths and calorimetric averaged",100.,7.,9.);
	for(int i=0;i<1000;i++) {
		hdoublecandidate[i] = new TH1F(Form("hdoublecandidate%d",i),"Double peak candidate",timetemp,0,timetemp);
	}
	for(int i=0;i<20;i++) {
		hSilicon[i]= new TH1F(Form("hSilicon%d",i),"Silicon signal",4096,0,4095);
		hSilicong[i]= new TH1F(Form("hSilicong%d",i),"Silicon signal - frontback mult=1",4096,0,4095);
		hSilicongg[i]= new TH1F(Form("hSilicongg%d",i),"Silicon signal - frontback mult=1, no MM",500,0,10.);
	}
	hbeamEl = new TH2F("hbeamEl","12N half event energy vs length",100,0,1000000,64,0,320);
	hpeaks = new TH1F("hpeaks","Number of peaks in waveform",8,0,8);
	hbeamtime = new TH2F("hbeamtime","Beam energy vs D2P time (check for register decay)",256,0,16384,100,0.,3e7);
	hbeamstop = new TH2F("hbeamstop","Beam stopping point",140,0,140,128,0.,128);
	hrndwalk = new TH1F("hrndwalk","Distance (xy) from end of beam to the centre of the decay",256,0,256);
	hrndwalkx = new TH1F("hrndwalkx","Distance (x) from end of beam to the centre of the decay",256,0,256);
	hrndwalky = new TH1F("hrndwalky","Distance (y) from end of beam to the centre of the decay",256,0,256);
	hrndwalktime = new TH2F("hrndwalktime","Distance (xy) from end of beam to the centre of the decay vs D2P time (ms)",100,0,50,100,0,50);
	hrndwalk2D = new TH2F("hrndwalk2D","Distance (xy) from end of beam to the centre of the decay x vs y",100,0,50,100,0,50);
	htheta2 = new TH1F("htheta2","Angle of 2nd arm relative to long arm",90,0,180);
	htheta2kf = new TH1F("htheta2kf","Angle of 2nd arm relative to long arm",90,0,180);
	htheta2theta3 = new TH2F("htheta2theta3","Angle of 2nd and 3rd arm (degrees)",45,0,180,45,0,180);
	htheta2theta3fit = new TH2F("htheta2theta3fit","Angle of 2nd and 3rd arm (degrees)",45,0,180,45,0,180);
	htheta2theta3kf = new TH2F("htheta2theta3kf","Angle of 2nd and 3rd arm (degrees) kinematically fitted",45,0,180,45,0,180);
	htheta2theta3theory = new TH2F("htheta2theta3theory","Angle of 2nd and 3rd arm (degrees) theory",45,0,180,45,0,180);
	hDalitz = new TH2F("hDalitz","Dalitz plot",100,-.5,.5,100,-.5,.5);
	hDalitztheory = new TH2F("hDalitztheory","Dalitz plot",100,-.5,.5,100,-.5,.5);
	hDalitzkf = new TH2F("hDalitzkf","Dalitz plot kinematically fitted",100,-1.,1.,100,-1.,1.);
	hthetapsi = new TH2F("hthetapsi","Direction of long arm vector",36,0,180,36,-180,180);
	hbeamEdE = new TH2F("hbeamEdE","Energy is front pad against the distance the beam travelled",100,0,10000,128,0,128);
	hEtotntot = new TH2F("hEtotntot","Total energy against total hits",100,0,1e7,100,0,200);
	hD2Ptime = new TH1F("hD2Ptime","D2P gated on Hoyle",100,0,100);
	harm_height = new TH2F("harm_height","Arm length versus height travelled",50,0.,150,50,0.,100.);
	hHoughXYZ = new TH3F("hHoughXYZ","3D Hough",180,-100,100,180,0,180,180,0,180);
	hfitchisqrd = new TH1F("hfitchisqrd","Average distance to point",100,0,4);
	hbdivyz = new TH1F("hbdivyz","Divergence of beam in YZ plane (theta)",100,-45,45);
	hbdivxy = new TH1F("hbdivxy","Divergence of beam in XY plane (theta)",100,-45,45);
	hdiv = new TH2F("hdiv","Divergence of decayno 0 vs decayno1 (theta vs theta)",100,0,180,100,0,180);
	hallsihits = new TH2F("hallsihits","All silicon hits",10,0,10,4,0,4);
	hDSSD = new TH2F("hDSSD","DSSD hits",16,0,16,16,0,16);
	hdriftcheck = new TH2F("hdriftcheck","Drift check (xy) vs (yz)",50,0,100,50,0,100);
	hwidth = new TH1F("hwidth","Width of deconvolved waveform",100,0,10);
	hdist12 = new TH2F("hdist12","Distance to arm 1 vs distance to arm 2",50,0,100,50,0,100);
	hrangeE1 = new TH2F("hrangeE1","Range of 1st alpha vs total signal",50,0,100,50,0,20000);
	hrangeE2 = new TH2F("hrangeE2","Range of 2nd alpha vs total signal",50,0,100,50,0,20000);
	hrangeE3 = new TH2F("hrangeE3","Range of 3rd alpha vs total signal",50,0,100,50,0,20000);
	hrange12 = new TH2F("hrange12","Range of 1st alpha vs 2nd alpha",50,0,100,50,0,100);
	hrange23 = new TH2F("hrange23","Range of 2nd alpha vs 3rd alpha",50,0,100,50,0,100);
	hExl = new TH2F("hExl","Excitation vs range of long arm",50,7,9,80,0,80);
	hEtotal = new TH1F("hEtotal","Total energy of 3 alphas from range",100,0,0.5);
	hpx_kf = new TH1F("hpx_kf","Total momentum in X",100,-0.5,0.5);
	hpy_kf = new TH1F("hpy_kf","Total momentum in Y",100,-0.5,0.5);
	hpz_kf = new TH1F("hpz_kf","Total momentum in Z",100,-0.5,0.5);
	hpx = new TH1F("hpx","Total momentum in X",100,-0.5,0.5);
	hpy = new TH1F("hpy","Total momentum in Y",100,-0.5,0.5);
	hpz = new TH1F("hpz","Total momentum in Z",100,-0.5,0.5);
	hEventID = new TH1F("hEventID","Event identification function",500,-1.,1.);
	hEventIDgoodevent = new TH1F("hEventIDgoodevent","Event identification function for good/bad event",500,0.,50.);
	hEventIDcut = new TH1F("hEventIDcut","Event identification function chi-sqrd<10",500,-1.,1.);
	hEventID2D = new TH2F("hEventID2D","Event identification function 2D",100,0,100.,100,0.,100);
	hEventIDheight = new TH2F("hEventIDheight","Event identification function vs height spread",100,-1,1,100,0.,200);
	hEventIDlength = new TH2F("hEventIDlength","Event identification function vs length of long arm",100,-1,1,100,0.,50);
	hEventIDposx = new TH2F("hEventIDposx","Event identification function vs maxx",100,-1,1,140,0.,140);
	hWiggle1 = new TH2F("hWiggle1","Range of alpha vs Wiggliness",50,0,50,50,0,15);
	hWiggle2 = new TH2F("hWiggle2","Range of alpha vs Wiggliness",50,0,50,50,0,15);
	hWiggle3 = new TH2F("hWiggle3","Range of alpha vs Wiggliness",50,0,50,50,0,15);
	hWiggleL1 = new TH2F("hWiggleL1","Displacement vs Wiggliness",50,0,50,50,0,15);
	hWiggleL2 = new TH2F("hWiggleL2","Displacement vs Wiggliness",50,0,50,50,0,15);
	hWiggleL3 = new TH2F("hWiggleL3","Displacement vs Wiggliness",50,0,50,50,0,15);
	hcurve = new TH2F("hcurve","Old length of arm 1 vs new length of curved arm",100,0,50,100,0,50);
	hcolinear = new TH1F("hcolinear","Colinearity of 3 alpha decay",100,-1,1);
	hbeamstopx = new TH1F("hbeamstopx","Beam stop x",140,0,140);
	hbeamstopy = new TH1F("hbeamstopy","Beam stop y",128,0,128);
	hdisplacement = new TH1F("hdisplacement","Displacement between beam stop and FITTED center",100,0,20);
	hRtotRmax = new TH2F("hRtotRmax","Total arm against longest arm",150,0,150,100,0,50);
	hHoylestopY = new TH2F("hHoylestopY","Hoyle beam stop",140,0,140,128,0,128);
	hd2ptimebeamy = new TH2F("hd2ptimebeamy","d2ptime vs beam y stop",30,0,30,128,0,128);
	hEtheta = new TH2F("hEtheta","C12 E vs theta",90,0,90,60,0,60);
	hEthetatheory = new TH2F("hEthetatheory","C12 E vs theta theory",90,0,90,100,0,5);
	hneutronenergy = new TH1F("hneutronenergy","Reconstructed neutron energy",1500,0,15);
	hneutronenergy2D = new TH2F("hneutronenergy2D","Reconstructed neutron energy",1500,0,15,100,0,40);
	hscatterEl = new TH2F("hscatterEl","Energy vs length",200,0,1e6,200,0,50);
	hErEs = new TH2F("hErEs","Energy sum of alpha from range vs integrated signal",100,0,2,100,0,5e6);
	hE1E2 = new TH2F("hE1E2","Energy integral vs peak",100,0,5e6,100,0,1e6);
	hlinearity = new TH2F("hlinearity","Linearity (first Hough vs second)",100,0,20,100,0,10);
	hChi1 = new TH1F("hChi1","Chi linear",100,0,10);
	hChi2 = new TH1F("hChi2","Chi 3alpha",100,0,10);
	hQvalue12C = new TH1F("hQvalue12C","Q-value for (n,a) assuming 13C and alpha",100,0,10);
	hQvalue16O = new TH1F("hQvalue16O","Q-value for (n,a) assuming 13C and alpha",100,0,10);
	hQvalue2D = new TH2F("hQvalue2D","Q-value for (n,a)",150,0,15,150,0,15);
	hthetatheta = new TH2F("hthetatheta","Theta long vs theta short",90,0,180,90,0,180);
	hthetathetaraw = new TH2F("hthetathetaraw","Theta long vs theta short",90,0,180,90,0,180);
	hthetatheta12C = new TH2F("hthetatheta12C","Theta long vs theta short",90,0,180,90,0,180);
	hthetatheta16O = new TH2F("hthetatheta16O","Theta long vs theta short",90,0,180,90,0,180);
	hthetatheta16Op = new TH2F("hthetatheta16Op","Theta long vs theta short",90,0,180,90,0,180);
	hthetathetatheory = new TH2F("hthetathetatheory","Theta long vs theta short theory",90,0,180,90,0,180);
	hthetathetatheoryCM = new TH2F("hthetathetatheoryCM","Theta lab vs theta CM theory",90,0,180,90,0,180);
	hEavg = new TH2F("hEavg","average E/length for long short",100,0,4000,100,0,4000);
	hchisqrdnalpha = new TH1F("hchisqrdnalpha","",500,0,100);
	hchi12Ctheta = new TH2F("hchi12Ctheta","",100,0,50,90,0,180);
	hchi16Otheta = new TH2F("hchi16Otheta","",100,0,50,90,0,180);
	hchi16Optheta = new TH2F("hchi16Optheta","",100,0,50,90,0,180);
	helscat = new TH2F("helscat","",180,-90,90,400,0,400);
	hlightR =new TH2F("hlightR","",100,0,5,100,0,5);
	hheavyR =new TH2F("hheavyR","",100,0,5,100,0,5);
//	hwidth = new TH1F("hwidth","Width of largest Gaussian",100,0,50);
//	hHoughPlane1 = new TH3F("hHoughPlane1","3D Hough normal psi vs theta vs x'",3600,0,360,900,-45,45,140,-122.5,122.5);
//	hHoughPlane2 = new TH3F("hHoughPlane2","3D Hough normal psi vs theta vs y'",3600,0,360,900,-45,45,140,-122.5,122.5);
}
Analysis::~Analysis() {
	delete h1;
	delete h2;
	for(int i=0;i<c_events;i++) delete htrack3D[i];
	delete hcurrenttrack3D;
}
void Analysis::DefineHists() {

	return;
}

int Analysis::ExecuteList(TString filename) {
	ifstream inlist;
	globaleventno=0;
        inlist.open(Form(filename));
	if(!inlist) {
		cout<<BOLDRED<<"List of files: "<<filename<<" not opened!\n"<<RESET;
		return -1;
	}
	else cout<<BOLDMAGENTA<<"List of files: "<< filename<<" opened\n"<<RESET;
	string newfile,lastfile="";
        while (inlist.good()){
		inlist>>newfile;
		if(newfile!=lastfile){
			Execute(newfile.c_str());//check to see if the last file is run twice for some reason
			cout<<GREEN<<"DONE"<<RESET<<endl;
		}
		lastfile=newfile;
	}
	return 0;
}

void Analysis::ReadMap() {

//////////////////Read MM channel map////////////////////////////////////////////////////////////////////////////////////////////////
	ifstream inmap;
        inmap.open(Form("mapchantomm.txt"));
        int numtel,good;
	if(!inmap) cout<<BOLDRED<<"MM Mapping file not opened!\n"<<RESET;
	else cout<<GREEN<<"MM mapping file opened\n"<<RESET;
	int incobo,inasad,inaget,inchan,outpixx,outpixy;

        while (inmap.good()){
                inmap>>inasad>>inaget>>inchan>>outpixx>>outpixy;
		map_x[inasad][inaget][inchan]=outpixx;
		map_y[inasad][inaget][inchan]=outpixy;
        }
	return;
//////////////////Read si channel map////////////////////////////////////////////////////////////////////////////////////////////////
	ifstream inmap2;
        inmap2.open(Form("mapchantosi.txt"));
	if(!inmap2) cout<<BOLDRED<<"Si mapping file not opened!\n"<<RESET;
	else cout<<GREEN<<"Si mapping file opened\n"<<RESET;
        while (inmap2.good()){
                inmap2>>inasad>>inaget>>inchan>>outpixx>>outpixy;
		simap_x[inasad][inaget][inchan]=outpixx;
		simap_y[inasad][inaget][inchan]=outpixy;
        }
	return;
}

void Analysis::SetResponse() {
	TFile *f = new TFile("waveform_ohio.root");
	TCanvas *c1 = (TCanvas*) f->Get("c_9fc3630_projection_102");
	TH1D *hresponse =  (TH1D*) c1->GetPrimitive("slice_py_of_hwaveform1");
	for(int i=0;i<256;i++) {
		response[i]=double(hresponse->GetBinContent(i+350)-hresponse->GetBinContent(hresponse->GetMinimum()));
	}
	for(int i=256;i<512;i++) {
		response[i]=0.;

	}
	return;
}
void Analysis::SetResponse2() {
	ifstream inlist;
        inlist.open("pulseshape.txt");
	if(!inlist) {
		cout<<BOLDRED<<"Response file not opened!\n"<<RESET;
	}
	else cout<<BOLDMAGENTA<<"Response file opened\n"<<RESET;
	for(int t=0;t<timetemp;t++) response[t]=0;
	int t,resp,timed=0;
        while (inlist.good()){
		inlist>>resp;
		response[timed]=resp;
//		cout<<"RAW: "<<timed<<"\t"<<response[timed]<<endl;
		timed++;
	}
	return;
}
void Analysis::SetResponse3() {
	TFile *f = new TFile("pulser_ohio.root");
	TCanvas *c1 = (TCanvas*) f->Get("c_56540d7bbb00_projection_40002");
//	c1->ls();
	TH1D *hresponse =  (TH1D*) c1->GetPrimitive("slice_py_of_hwaveform360");
	for(int i=0;i<70;i++) {
		response[i]=double(hresponse->GetBinContent(i+440)-18500)*0.02;
	}
	for(int i=70;i<512;i++) {
		response[i]=0.02*2900.*TMath::Exp(-0.1*(i-70.));//phase out

	}
//	for(int i=0;i<512;i++) cout<<i<<"\t"<<response[i]<<endl;
	return;
}

void Analysis::HoughTransform(int mode) {//0 for beam only
	double theta=0.,costheta=0.,sintheta=0.,xdist=0.,ydist=0.,zdist=0.,radiusxy=0.,radiusyz=0.,radiusxz=0.;
	int xpix=0,ypix=0,zpix=0;
//	double scaletime=driftv*100;//1.1mm/timebucket
	double x_center =centr_dist_x;
	double y_center=centr_dist_y;
	double z_center = centr_dist_z;
//	double x_center =0.;
//	double y_center=0.;
//	double z_center = 70.;
	Dist(dmaxx,dmaxy,dmaxz,x_center,y_center,z_center);
	mmHit=128;
	TRandom3 *rndmgen = new TRandom3();
	int goodhits=0;
//	hcurrenttrack3D->Fill(x_center,y_center,z_center,1e4);
	for(int theta_i=0;theta_i<1800;theta_i++) {
		theta=0.1*theta_i;
		costheta=cos(theta*d2r);
		sintheta=sin(theta*d2r);
		double xdist,ydist,zdist;
		xdist=rnd3->Rndm()-0.5;//-0.5->0.5
		ydist=rnd3->Rndm()-0.5;
		zdist=rnd3->Rndm()-0.5;
		radiusxy=xdist*costheta+ydist*sintheta;
		radiusyz=ydist*costheta+zdist*sintheta;
		radiusxz=xdist*costheta+zdist*sintheta;
		hHoughXYb->Fill(theta,radiusxy,1e5);
		hHoughYZb->Fill(theta,radiusyz,1e5);
		hHoughXZb->Fill(theta,radiusxz,1e5);
						if(mode==2 && good_evt<draw_events) {
							hHoughXY[good_evt]->Fill(theta,radiusxy,1e5);
							hHoughYZ[good_evt]->Fill(theta,radiusyz,1e5);
							hHoughXZ[good_evt]->Fill(theta,radiusxz,1e5);
						}

	}

	for(int x=0;x<140;x++) {
		for(int y=0;y<128;y++) {
			for(int z=0;z<timetemp;z++) {
//				if((hcurrenttrack3D->GetBinContent(x,y,z)>0)) {//choose hbeamtrack for beam and hcurrenttrack3D for decay
				if((hcurrenttrack3D->GetBinContent(x+1,y+1,z+1)>0)) {//choose hbeamtrack for beam and hcurrenttrack3D for decay
					double maxE;
					goodhits++;
//					maxE=hcurrenttrack3D->GetBinContent(x,y,z);
					maxE=hcurrenttrack3D->GetBinContent(x+1,y+1,z+1);
					xpix=x;
					ypix=y;
					zpix=z;
					Dist(xpix,ypix,zpix,xdist,ydist,zdist);
					if(good_evt<draw_events) hdistxy[good_evt]->Fill(xdist,ydist);
					if(good_evt<draw_events) hdistxz[good_evt]->Fill(xdist,zdist);
					if(good_evt<draw_events) hdistyz[good_evt]->Fill(ydist,zdist);
//
					if(good_evt<draw_events) {
						if(mode==1) {
							hdistxy1[good_evt]->Fill(xdist,ydist);
							hdistxz1[good_evt]->Fill(xdist,zdist);
							hdistyz1[good_evt]->Fill(ydist,zdist);
						}
						if(mode==2) {
							hdistxy2a[good_evt]->Fill(xdist,ydist);
							hdistxz2a[good_evt]->Fill(xdist,zdist);
							hdistyz2a[good_evt]->Fill(ydist,zdist);
							hdistxy2b[good_evt]->Fill(xdist,ydist);
							hdistxz2b[good_evt]->Fill(xdist,zdist);
							hdistyz2b[good_evt]->Fill(ydist,zdist);
							hdistxy2c[good_evt]->Fill(xdist,ydist);
							hdistxz2c[good_evt]->Fill(xdist,zdist);
							hdistyz2c[good_evt]->Fill(ydist,zdist);
						}
						if(mode==3) {
							hdistxy3[good_evt]->Fill(xdist,ydist);
							hdistxz3[good_evt]->Fill(xdist,zdist);
							hdistyz3[good_evt]->Fill(ydist,zdist);
						}
					}
//
					xdist-=x_center;
					ydist-=y_center;
					zdist-=z_center;
					for(int theta_i=0;theta_i<1800;theta_i++) {
						theta=0.1*theta_i;
						costheta=cos(theta*d2r);
						sintheta=sin(theta*d2r);
						radiusxy=xdist*costheta+ydist*sintheta;
						radiusyz=ydist*costheta+zdist*sintheta;
						radiusxz=xdist*costheta+zdist*sintheta;
						hHoughXYb->Fill(theta,radiusxy,maxE);
						hHoughYZb->Fill(theta,radiusyz,maxE);
						hHoughXZb->Fill(theta,radiusxz,maxE);
						if(mode==2 && good_evt<draw_events) {
							hHoughXY[good_evt]->Fill(theta,radiusxy,maxE);
							hHoughYZ[good_evt]->Fill(theta,radiusyz,maxE);
							hHoughXZ[good_evt]->Fill(theta,radiusxz,maxE);
						}
					}
				}
			}
		}
	}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double stddevxy,thetaxy=0.,zbin,ybin,xbin,thetayz=0.,stddevyz,thetaxz=0.,stddevxz;
	double minstddevxy=100000,minstddevyz=100000,minstddevxz=100000;
	int entries=hHoughXYb->ProjectionY("projy",450,451)->GetEntries();//not sure why 450,451-- 45 degrees
	if(goodhits<3) {
		length[mode-1]=0.;
		palpha[mode-1][0]=0.;
		palpha[mode-1][1]=0.;
		palpha[mode-1][2]=0.;
		hHoughXYb->Reset();
		hHoughYZb->Reset();
		hHoughXZb->Reset();
		return;
	}
        double cur_max=0;
        for(int x=0;x<hHoughXYb->GetNbinsX();x++) {
                for(int y=0;y<hHoughXYb->GetNbinsY();y++) {
                        if(hHoughXYb->GetBinContent(x+1,y+1)>cur_max) {
//			cout<<x<<"\t"<<y<<"\t"<<cur_max<<endl;

                                cur_max=hHoughXYb->GetBinContent(x+1,y+1);
                                thetaxy=hHoughXYb->GetXaxis()->GetBinCenter(x);
                                radiusxy=hHoughXYb->GetYaxis()->GetBinCenter(y);
                        }
                }
        }
        cur_max=0;
        for(int x=0;x<hHoughYZb->GetNbinsX();x++) {
                for(int y=0;y<hHoughYZb->GetNbinsY();y++) {
                        if(hHoughYZb->GetBinContent(x+1,y+1)>cur_max) {
                                cur_max=hHoughYZb->GetBinContent(x+1,y+1);
                                thetayz=hHoughYZb->GetXaxis()->GetBinCenter(x);
                                radiusyz=hHoughYZb->GetYaxis()->GetBinCenter(y);
                        }
                }
        }
        cur_max=0;
        for(int x=0;x<hHoughXZb->GetNbinsX();x++) {
                for(int y=0;y<hHoughXZb->GetNbinsY();y++) {
                        if(hHoughXZb->GetBinContent(x+1,y+1)>cur_max) {
                                cur_max=hHoughXZb->GetBinContent(x+1,y+1);
                                thetaxz=hHoughXZb->GetXaxis()->GetBinCenter(x);
                                radiusxz=hHoughXZb->GetYaxis()->GetBinCenter(y);
                        }
                }
        }
//	cout<<"XY: "<<thetaxy<<"\t"<<radiusxy<<endl;
////////////////
	double costhetaxy=cos(thetaxy*d2r);
	double sinthetaxy=sin(thetaxy*d2r);
	double ifzeroxy=sinthetaxy*costhetaxy;
	if(ifzeroxy==0) {
		theta+=0.1;
		costhetaxy=cos(thetaxy*d2r);
		sinthetaxy=sin(thetaxy*d2r);
	}
	double maxbinxz = hHoughXZb->GetMaximumBin();
	double costhetaxz = cos(thetaxz*d2r);
	double sinthetaxz = sin(thetaxz*d2r);
	double ifzeroxz=sinthetaxz*costhetaxz;
	if(ifzeroxz==0) {
                theta+=0.1;
                costhetaxz=cos(thetaxz*d2r);
                sinthetaxz=sin(thetaxz*d2r);
        }
	double maxbinyz = hHoughYZb->GetMaximumBin();
        double costhetayz = cos(thetayz*d2r);
        double sinthetayz = sin(thetayz*d2r);
        double ifzeroyz=sinthetayz*costhetayz;
        if(ifzeroyz==0) {
                theta+=0.1;
                costhetayz=cos(thetayz*d2r);
                sinthetayz=sin(thetayz*d2r);
        }
///////////////////////////////////////
	double posx,posy,posz,timez,timemm;
	posx=(radiusxy+300*sinthetaxy)/costhetaxy;
	double xline[460],yline[460],zline[460];
	int loopc=0;
	for(int j=-224;j<224;j++) {
		posx=(radiusxy-j*sinthetaxy)/costhetaxy;
		timemm=((radiusyz-j*costhetayz)/sinthetayz)*100.;//100ns per time bin
		posz= (timemm*driftv);
		posz= (((radiusyz-j*costhetayz)/sinthetayz));
		if(good_evt<draw_events) hdistxy[good_evt]->Fill(posx+x_center,j+y_center,1e5*mode);
		if(good_evt<draw_events) hdistyz[good_evt]->Fill(j+y_center,posz+z_center,1e5*mode);
		if(good_evt<draw_events) hdistxz[good_evt]->Fill(posx+x_center,posz+z_center,1e5*mode);
		if(good_evt<draw_events) {
			if(mode==1) {
				if(good_evt<draw_events) hdistxy1[good_evt]->Fill(posx+x_center,j+y_center,1e5*mode);
				if(good_evt<draw_events) hdistyz1[good_evt]->Fill(j+y_center,posz+z_center,1e5*mode);
				if(good_evt<draw_events) hdistxz1[good_evt]->Fill(posx+x_center,posz+z_center,1e5*mode);

			}
			if(mode==2) {
				if(good_evt<draw_events) hdistxy2a[good_evt]->Fill(posx+x_center,j+y_center,1e5*mode);
				if(good_evt<draw_events) hdistyz2a[good_evt]->Fill(j+y_center,posz+z_center,1e5*mode);
				if(good_evt<draw_events) hdistxz2a[good_evt]->Fill(posx+x_center,posz+z_center,1e5*mode);

			}
			if(mode==3) {
				if(good_evt<draw_events) hdistxy3[good_evt]->Fill(posx+x_center,j+y_center,1e5*mode);
				if(good_evt<draw_events) hdistyz3[good_evt]->Fill(j+y_center,posz+z_center,1e5*mode);
				if(good_evt<draw_events) hdistxz3[good_evt]->Fill(posx+x_center,posz+z_center,1e5*mode);

			}
		}
		xline[loopc]=posx+x_center;
		yline[loopc]=j+y_center;
		zline[loopc]=posz+z_center;
//		cout<<xline[loopc]<<"\t"<<yline[loopc]<<"\t"<<zline[loopc]<<"\t"<<loopc<<endl;
		loopc++;
	}
	int maxj=0;
	double thislength=0.,maxlength=0.;
	double endat_x=0,endat_y=0,endat_z=0,totaldx=0.,totaldy=0.,totaldz=0.;
//	cout<<centr_dist_y<<"\t"<<dmaxy<<endl;
	int removed=0;
	remove[mode]=0;

	for(int x=0;x<140;x++) {
		for(int y=0;y<128;y++) {
			for(int z=0;z<timetemp;z++) {
				if((hcurrenttrack3D->GetBinContent(x+1,y+1,z+1)>0)) {//choose hbeamtrack for beam and hcurrenttrack3D for decay
					for(int j=0;j<2*224+1;j++) {//loop over points on y line
						double distance,atx,aty,atz;
						Dist(x,y,z,atx,aty,atz);
						distance=sqrt(TMath::Power(atx-xline[j],2.)+TMath::Power(aty-yline[j],2.)+TMath::Power(atz-zline[j],2.));
//						cout<<xline[j]<<"\t"<<yline[j]<<"\t"<<zline[j]<<"\t"<<atx<<"\t"<<aty<<"\t"<<atz<<"\t"<<distance<<"\t"<<((1.+0.1*fabs(aty-centr_dist_y)))<<endl;

//						if(distance<((1.+0.1*fabs(aty-centr_dist_y))) && (fexperiment==1 || ((mode==1 && j>224) || mode==2))) {//find the forward arm first
						if(distance<5. && (fexperiment==1 || ((mode==1 && j>224) || mode==2))) {//find the forward arm first
//							cout<<"removed"<<endl;
							removed++;
							remove[mode-1]++;
							hcurrenttrack3D->SetBinContent(x+1,y+1,z+1,-fabs(hcurrenttrack3D->GetBinContent(x+1,y+1,z+1)));//remove this as close to current line
							thislength=sqrt(TMath::Power(atx-x_center,2.)+TMath::Power(aty-y_center,2.)+TMath::Power(atz-z_center,2.));//length
							if(thislength>maxlength) {
								maxlength=thislength;
								endat_x=x;
								endat_y=y;
								endat_z=z;
								totaldx=atx-x_center;
								totaldy=aty-y_center;
								totaldz=atz-z_center;
							}
						}
					}
				}
			}
		}
	}
	if(mode==2) hlinearity->Fill(remove[0],remove[1]);


	if(good_evt<draw_events) {
		if(mode==1) hdistxy1[good_evt]->SetTitle(Form("Removed: %d",removed));
		if(mode==2) hdistxy2a[good_evt]->SetTitle(Form("Removed: %d",removed));
		if(mode==3) hdistxy3[good_evt]->SetTitle(Form("Removed: %d",removed));
	}
//	cout<<maxlength<<endl;
	length[mode-1]=maxlength;
	palpha[mode-1][0]=totaldx;
	palpha[mode-1][1]=totaldy;
	palpha[mode-1][2]=totaldz;
	scatterangle=atan(sqrt(totaldx*totaldx+totaldz*totaldz)/totaldy);
	double divxy,divyz;
	divxy=thetaxy;
	if(divxy>90) divxy=180-divxy;//0->90
	divyz=thetayz-90;
	if(divyz>90) divyz=180-divyz;//0->90
	if(good_evt<draw_events) htrackxy[good_evt]->Fill(endat_x,endat_y,1e6);
	if(good_evt<draw_events) htrackyz[good_evt]->Fill(endat_y,endat_z,1e6);
	if(good_evt<draw_events) htrackxz[good_evt]->Fill(endat_x,endat_z,1e6);
	hHoughXYb->Reset();
	hHoughYZb->Reset();
	hHoughXZb->Reset();
	return;
}


void Analysis::SetExperiment(int setvar) {

	fexperiment=setvar;
	return;
}


int Analysis::HoughTransformEvent() {
        double theta=0.,costheta=0.,sintheta=0.,xdist=0.,ydist=0.,zdist=0.,radiusxy=0.,radiusyz=0.,radiusxz=0.,totalradius=0.;
        double d2r=atan(1.)/45.;//4*atan(1.)=pi-->pi/180=atan(1.)/45
        int xpix=0,ypix=0,zpix=0;
//        double scaletime=1.1;//1.1mm/timebucket
        mmHit=128;
        TRandom3 *rndmgen = new TRandom3();
	int counter=0;
	double a1,b1,E;
	a1=rndmgen->Rndm();
	b1=rndmgen->Rndm();
	double x_center,y_center,z_center;
	Dist(decay_x,decay_y,decay_z,x_center,y_center,z_center);//get the distances for the center of the decay and shift everything relative to this
	const int smeart=100;
	double center_weighting=2000;//How "hard" do you want to drag the lines towards the decay center given? 2000 usually forces everything through it
	bool get_lines=true;//set to true until the lines are a good enough fit
	double total_three=0.;
	const int numb_weightings=5;
	double rxy[3][numb_weightings],ryz[3][numb_weightings],rxz[3][numb_weightings],rthetaxy[3][numb_weightings],rthetayz[3][numb_weightings],rthetaxz[3][numb_weightings],chisqrd[numb_weightings];
	int goodhits=0;
	y_center=centr_dist_y;
	z_center=centr_dist_z;
	hHoughYZ2 = new TH2F("hHoughYZ2","",1800,0,180,600,-50,50);
	for(int y=0;y<128;y++) {
		for(int z=0;z<timetemp;z++) {
			if((hcurrentYZ->GetBinContent(y+1,z+1)>0)) {//choose hbeamtrack for beam and hcurrenttrack3D for decay
				double maxE;
				goodhits++;
				maxE=hcurrentYZ->GetBinContent(y+1,z+1);
				ypix=y;
				zpix=z;
				Dist(0,ypix,zpix,xdist,ydist,zdist);
				ydist-=y_center;
				zdist-=z_center;
//				cout<<ydist<<"\t"<<zdist<<endl;
				for(int theta_i=0;theta_i<1800;theta_i++) {
					theta=0.1*theta_i;
					costheta=cos(theta*d2r);
					sintheta=sin(theta*d2r);
					radiusyz=ydist*costheta+zdist*sintheta;
//					hHoughYZ2->Fill(theta,radiusyz,maxE);
					hHoughYZ2->Fill(theta,radiusyz);
					if(good_evt<draw_events) {
//						hdistyz[good_evt]->Fill(ydist+y_center,zdist);
						hHoughXY2[good_evt]->Fill(theta,radiusyz,maxE);
					}
				}
			}
		}
	}
	double maxXYZ=0.,max_X=0,max_Y=0,max_Z=0;
////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double stddevxy,thetaxy=0.,zbin,ybin,xbin,thetayz=0.,stddevyz,thetaxz=0.,stddevxz;
        double minstddevxy=100000,minstddevyz=100000,minstddevxz=100000;
	int binxy=0,binyz=0;
	if(goodhits<3) {
		delete hHoughYZ2;
		return 0;
	}
	double xline[460],yline[460],zline[460];

        double cur_max=0;
        for(int x=0;x<hHoughYZ2->GetNbinsX();x++) {
//                for(int y=0.5*hHoughYZ2->GetNbinsY()-1;y<0.5*hHoughYZ2->GetNbinsY()+1;y++) {//only search around vertex
                for(int y=0;y<hHoughYZ2->GetNbinsY();y++) {//only search around vertex
                        if(hHoughYZ2->GetBinContent(x,y)>cur_max) {
                                cur_max=hHoughYZ2->GetBinContent(x,y);
                                thetayz=0.5*(hHoughYZ2->GetXaxis()->GetBinCenter(x)+hHoughYZ2->GetXaxis()->GetBinCenter(x+1));
                                radiusyz=0.5*(hHoughYZ2->GetYaxis()->GetBinCenter(y)+hHoughYZ2->GetYaxis()->GetBinCenter(y+1));
                        }
                }
        }
//	radiusyz=0.;//Force to go to vertex
	cout<<"theta /radius "<<thetayz<<"\t"<<radiusyz<<endl;
	double costhetayz=cos(d2r*thetayz);
	double sinthetayz=sin(d2r*thetayz);
////////////////
	double posz,timemm;
	int loopc=0;
        for(int j=-224;j<224;j++) {
//                timemm=((radiusyz-j*costhetayz)/sinthetayz)*100.;//100ns per time bin
  //              posz= (timemm*driftv);
                posz= (((radiusyz-j*costhetayz)/sinthetayz));
                yline[loopc]=(j)+y_center;
                zline[loopc]=posz;
//		if(good_evt<draw_events) {
//			hdistyz[good_evt]->Fill(yline[loopc],zline[loopc]);
//		}
                loopc++;
        }

////////////////
		for(int y=0;y<128;y++) {
			for(int z=0;z<timetemp;z++) {
				if((hcurrentYZ->GetBinContent(y+1,z+1)>0)) {//choose hbeamtrack for beam and hcurrenttrack3D for decay
					for(int j=0;j<2*224+1;j++) {//loop over points on y line
						double distance,atx,aty,atz;
						Dist(0,y,z,atx,aty,atz);
						distance=sqrt(TMath::Power(aty-yline[j],2.)+TMath::Power(atz-zline[j],2.));
						if(distance<3.) {//
							hcurrentYZ->SetBinContent(y+1,z+1,0.);//remove this as close to current line
						}
					}
				}
			}
		}
	delete hHoughYZ2;

        return 1;
}


void Analysis::FindPlane() {
//	Get the reaction plane by generating a direction of an icosahedron and then do the dot product of all vectors in the point cloud with this vector and minimize
//	Generate the directions
//	12 vertices = {(0,pm1,pmr),(pm1,pmr,0),(pmr,0,pm1)} r = 0.5(1+sqrt(5))
	double n[3]={0.},theta,psi,bestdot=10000;
	int best_i=-1,best_j=-1;
	double best_theta=0.,best_psi=0.;
	const int theta_points=10,psi_points=10;
	double decay_dist_x,decay_dist_y,decay_dist_z;
	Dist(decay_x,decay_y,decay_z,decay_dist_x,decay_dist_y,decay_dist_z);
	for(int i=0;i<=theta_points;i++) {
//		theta=pi*i/(1.*theta_points);
		theta=acos((2.*i)/(1.*theta_points)-1.);
		for(int j=0;j<=psi_points;j++) {
			psi=2.*pi*j/(1.*psi_points);
			n[0]=sin(theta)*cos(psi);
			n[1]=sin(theta)*sin(psi);
			n[2]=cos(theta);
//			Loop over all data points
			double thisdot=0.;
			for(int x=0;x<140;x++) {
				for(int y=0;y<128;y++) {
					for(int z=0;z<timetemp;z++) {
						if(hcurrenttrack3D->GetBinContent(x+1,y+1,z+1)>0) {
//							Calculate dot product with the normal vector
							double xd=0.,yd=0.,zd=0.;
							Dist(x,y,z,xd,yd,zd);
							double p[3]={0.},pt;
							p[0]=xd-decay_dist_x;
							p[1]=yd-decay_dist_y;
							p[2]=zd-decay_dist_z;
							pt=sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
							for(int m=0;m<3;m++) p[m]/=pt;//Normalize
							thisdot+=(p[0]*n[0])+(p[1]*n[1])+(p[2]*n[2]);
						}//check for hit
					}//z
				}//y
			}//x
			if(fabs(thisdot)<bestdot) {//update best
				bestdot=fabs(thisdot);
				best_i=i;
				best_j=j;
				best_theta=theta;
				best_psi=psi;
			}
		}
	}
//	cout<<"BEST: "<<best_theta/d2r<<"\t"<<best_psi/d2r<<"        "<<sin(best_theta)*cos(best_psi)<<","<<sin(best_theta)*sin(best_psi)<<","<<cos(best_theta)<<endl;
	Normal[0]=sin(best_theta)*cos(best_psi);
	Normal[1]=sin(best_theta)*sin(best_psi);
	Normal[2]=cos(best_theta);
	return;
}


void Analysis::HoughTransform2() {
	double theta=0.,costheta=0.,sintheta=0.,xdist=0.,ydist=0.,zdist=0.,radiusxy=0.,radiusyz=0.,radiusxz=0.;
	int xpix=0,ypix=0,zpix=0;
//	double scaletime=1.1;//1.1mm/timebucket
	double x_center =122.5;
	TRandom3 *rndmgen = new TRandom3();
	TH2F *hHoughXYflat = new TH2F("hHoughXYflat","Hough on flattened decay",180,0,180,200,-100,100);
	for(int x=0;x<100;x++) {
		for(int y=0;y<100;y++) {
			if(hflat_all->GetBinContent(x+1,y+1)>0) {
				double maxE=0.,maxT=0.;
				double xdist,ydist;
				xdist=hflat_all->GetXaxis()->GetBinCenter(x);
				ydist=hflat_all->GetYaxis()->GetBinCenter(y);
				for(int theta_i=0;theta_i<1800;theta_i++) {
					theta=0.1*theta_i;
					costheta=cos(theta*d2r);
					sintheta=sin(theta*d2r);
					radiusxy=xdist*costheta+ydist*sintheta;
					hHoughXYflat->Fill(theta,radiusxy);
//					cout<<theta<<"\t"<<radiusxy<<endl;
				}
			}
		}
	}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double stddevxy,thetaxy=0.,zbin,ybin,xbin,thetayz=0.,stddevyz,thetaxz=0.,stddevxz;
	double minstddevxy=100000,minstddevyz=100000,minstddevxz=100000;
	int entries=hHoughXYflat->ProjectionY("projy",450,451)->GetEntries();//
/*
	for(int i=0;i<1800;i++) {
/////////////////////////////////////////////////////////////////////
///////////	Find sharpest projection point at 45 degrees	/////
/////////////////////////////////////////////////////////////////////
///	XY
		if(hHoughXYflat->ProjectionY("projy",i,i+1)->GetEntries()==entries) {
			stddevxy=hHoughXYflat->ProjectionY("projy",i,i+1)->GetStdDev();
//			cout<<stddevxy<<"\t"<<i<<endl;
			if(minstddevxy>=stddevxy) {
				minstddevxy=stddevxy;
				thetaxy=hHoughXYflat->GetXaxis()->GetBinCenter(i);
				ybin=hHoughXYflat->ProjectionY("projy",i,i+1)->GetMaximumBin();
				radiusxy=hHoughXYflat->GetYaxis()->GetBinCenter(ybin);
			}
		}
	}
*/
	double cur_max=0;
	for(int x=0;x<hHoughXYflat->GetNbinsX();x++) {
		for(int y=0;y<hHoughXYflat->GetNbinsY();y++) {
			if(hHoughXYflat->GetBinContent(x+1,y+1)>cur_max) {
				cur_max=hHoughXYflat->GetBinContent(x+1,y+1);
				thetaxy=hHoughXYflat->GetXaxis()->GetBinCenter(x);
                                radiusxy=hHoughXYflat->GetYaxis()->GetBinCenter(y);
			}
		}
	}
////////////////////////////////////////
	double costhetaxy=cos(thetaxy*d2r);
	double sinthetaxy=sin(thetaxy*d2r);
	double ifzeroxy=sinthetaxy*costhetaxy;
	if(ifzeroxy==0) {
		theta+=0.1;
		costhetaxy=cos(thetaxy*d2r);
		sinthetaxy=sin(thetaxy*d2r);
	}
//	cout<<evtno<<"\t"<<radiusxy<<"\t"<<thetaxy<<endl;
///////////////////////////////////////
	double posx,posy,posz,timez,timemm;
	for(int j=-300;j<400;j++) {
		posx=(radiusxy-j*sinthetaxy)/costhetaxy;
		if(good_evt<draw_events) hflat[good_evt]->Fill(posx,j);
	}
	delete hHoughXYflat;
	return;
}

void Analysis::SetHistName(TString fname) {
	cout<<"HIST NAME: "<<fname<<endl;
	outputdefaultname=fname;
	return;
}

void Analysis::HoughTransformGeneric(TH2F *hgeneric,double &radius, double &thetaout) {//feed a generic TH2F histogram
	double theta=0.,costheta=0.,sintheta=0.,xdist=0.,ydist=0.,zdist=0.,radiusxy=0.,radiusyz=0.,radiusxz=0.;
	int xpix=0,ypix=0,zpix=0;
	TRandom3 *rndmgen = new TRandom3();
	TH2F *hgenericHough = new TH2F("hgenericHough","Hough on generic",180,0,180,200,-100,100);
	for(int x=0;x<hgeneric->GetNbinsX();x++) {
		for(int y=0;y<hgeneric->GetNbinsY();y++) {
			if(hgeneric->GetBinContent(x+1,y+1)>0) {
				double maxE=0.,maxT=0.;
				double xdist,ydist;
				xdist=hgeneric->GetXaxis()->GetBinCenter(x);
				ydist=hgeneric->GetYaxis()->GetBinCenter(y);
				for(int theta_i=0;theta_i<1800;theta_i++) {
					theta=0.1*theta_i;
					costheta=cos(theta*d2r);
					sintheta=sin(theta*d2r);
					radiusxy=xdist*costheta+ydist*sintheta;
					hgenericHough->Fill(theta,radiusxy);
//					cout<<theta<<"\t"<<radiusxy<<endl;
				}
			}
		}
	}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double stddevxy,thetaxy=0.,zbin,ybin,xbin,thetayz=0.,stddevyz,thetaxz=0.,stddevxz;
	double minstddevxy=100000,minstddevyz=100000,minstddevxz=100000;
	int entries=hgenericHough->ProjectionY("projy",450,451)->GetEntries();//
	double cur_max=0;
//	Use the highest peak value to avoid "averaging" between multiple lines
	for(int x=0;x<hgenericHough->GetNbinsX();x++) {
		for(int y=0;y<hgenericHough->GetNbinsY();y++) {
			if(hgenericHough->GetBinContent(x+1,y+1)>cur_max) {
				cur_max=hgenericHough->GetBinContent(x+1,y+1);
				thetaxy=hgenericHough->GetXaxis()->GetBinCenter(x);
                                radiusxy=hgenericHough->GetYaxis()->GetBinCenter(y);
			}
		}
	}
////////////////////////////////////////
	double costhetaxy=cos(thetaxy*d2r);
	double sinthetaxy=sin(thetaxy*d2r);
	double ifzeroxy=sinthetaxy*costhetaxy;
	if(ifzeroxy==0) {
		theta+=0.1;
		costhetaxy=cos(thetaxy*d2r);
		sinthetaxy=sin(thetaxy*d2r);
	}
///////////////////////////////////////
	double posx,posy,posz,timez,timemm;
	for(int j=-300;j<400;j++) {
		posx=(radiusxy-j*sinthetaxy)/costhetaxy;
	}
	delete hgenericHough;
	radius=radiusxy;//feed to extraction variables
	thetaout=thetaxy;
	return;
}

void Analysis::FillEnergyRange() {
	ifstream infile;
	infile.open("SRIM20Torr.txt");
	int i=0;
	while(infile.good() && i<ERentries){
		double ekeV,dum1,dum2,range,dum3,dum4;
		infile>>ekeV>>dum1>>dum2>>range>>dum3>>dum4;
		EnergyRange[0][i]=ekeV;
		EnergyRange[1][i]=range;
		i++;
	}
	return;
}

void Analysis::SetRunNumber(int runn) {
	runnumber=runn;
//	pressure=GetPressure(FullRange[runn]);
	if(runn>=223 && runn<=225) Ebeam=7.2;//4.25 MeV - 100 Torr
	if(runn>=219 && runn<=222) Ebeam=7.55;//4.55 MeV - 100 Torr
	if(runn>=215 && runn<=218) Ebeam=7.75;//4.75 MeV - 100 Torr
	if(runn>=207 && runn<=214) Ebeam=8.15;//5.15 MeV - 100 Torr
	if(runn>=201 && runn<=205) Ebeam=8.36;//5.36 MeV - 100 Torr
	if(runn>=108 && runn<=111) Ebeam=8.36;//5.36 MeV
	if(runn>=228 && runn<=232) Ebeam=8.6;//5.6 MeV
	if(runn>=154 && runn<=170) Ebeam=8.75;//5.75
	if(runn>=233 && runn<=240) Ebeam=8.875;//5.875
	if((runn>=177 && runn<=182) || runn==267) Ebeam=9.;//6
	if(runn>=241 && runn<=249) Ebeam=9.125;//6.125
	if(runn>=87 && runn<=100) Ebeam=9.25;//6.25
	if(runn>=186 && runn<=193) Ebeam=9.5;//6.5
	if(runn>=140 && runn<=151) Ebeam=9.75;//6.75
	if(runn>=78 && runn<=85) Ebeam=10.;//7
//	pressure=20.;
	cout<<"Ebeam: "<<Ebeam<<" for run "<<runn<<endl;
	return;
}

double Analysis::pulseshape(Double_t *x, Double_t *p) {
        Double_t semiGaus = p[1]*22.68113723*TMath::Exp(-3.0*(x[0] - p[2])/p[3])*sin((x[0] - p[2])/p[3])*pow((x[0] - p[2])/p[3], p[4]);
        return (x[0] >= p[2]) ? p[0] + x[0]*p[5] + semiGaus : p[0] + x[0]*p[5];
}

void Analysis::FitSaturated(double& maxE,double& maxT,double centroid,TH1F* hfit) {
        TF1* f1 = new TF1("myfunc",this,&Analysis::pulseshape,0,timetemp,6);
        double width=8;
//	cout<<maxE<<"\t"<<maxT<<endl;
//	for(int i=0;i<512;i++) cout<<i<<"\t"<<hfit->GetBinContent(i+1)<<endl;
        f1->SetParName(0,"const BG");
        f1->SetParName(1,"Amp");
        f1->SetParName(2,"Centr");
        f1->SetParName(3,"Width");
        f1->SetParName(4,"Power");
        f1->SetParName(5,"linear BG");
        f1->SetParLimits(0,-200,200);
        f1->SetParLimits(1,10,10000);
        f1->SetParLimits(2,centroid-50-width,centroid+50+width);
        f1->SetParLimits(3,5.,200.);
        f1->SetParLimits(4,0.,5.);
        f1->SetParLimits(5,-1.,1.);
        f1->SetParameter(0,0.0000001);//bg
        f1->SetParameter(1,hfit->GetBinContent(centroid+1));//amplitude
        f1->SetParameter(2,centroid-width);//centroid
        f1->SetParameter(3,width);//width
        f1->SetParameter(4,1.45);//power
        f1->SetParameter(5,0.00000001);//bg linear
//      f1->Draw();
        hfit->Fit("myfunc","QNM");//N=do not draw, Q = quiet mode, M = More (improves fit when local minimum found)
//      hfit->Fit("myfunc");
        double chsqrd = f1->GetChisquare()/f1->GetNDF();//get the REDUCED chi-sqrd
  //      if(chsqrd>25.) {
//                if(verbose)cout<<BOLDRED<<"Fitting error on saturated waveform "<<RESET<<chsqrd<<endl;
/*
                TCanvas *c2 = new TCanvas("c2","",0,0,800,600);
                hfit->Draw();
                f1->Draw("SAME");
                c2->Update();
                for(double t=0;t<1000000000;t++) {
                        t=t+0.1;
                        if(chsqrd>10) t-=1.;
                }
                delete c2;
*/
  //	}
	for(int t=0;t<timetemp;t++) {
                hfit->SetBinContent(t+1,f1->Eval(t));//corrected waveform is the fitted waveform
        }
	double peakh=0;
//      peakh=f1->GetParameter(1);
	double peakl=0;
	if(chsqrd<25 && f1->GetMaximumX()>10) {
	        peakh=f1->GetMaximum();//gets the maximum of the fit using interpolation methods
		peakl=f1->GetMaximumX();
	}
	delete f1;
	maxE=peakh;
	maxT=peakl;
	return;
}

void Analysis::FillRunRange() {
        ifstream infile;
        infile.open("collected_lengths.txt");
        int i=0;
	for(int j=0;j<200;j++) {
		FullRange[j]=37.5;
	}
        while(infile.good() && i<ERentries){
		int run;
                double time,dum1,dum2,range,dum3,dum4,dum5;
                infile>>run>>time>>dum1>>dum2>>range>>dum3>>dum4>>dum5;
		FullRange[run]=range;
                i++;
        }
        return;
}

double Analysis::GetPressure(double range) {
	double pressure=50;

	return pressure;
}

TString Analysis::GetHistName() {
	cout<<"HIST NAME: "<<outputdefaultname<<endl;
	return outputdefaultname;
}


void Analysis::Ransac() {
//	Get point cloud
        vector<double> x;
        vector<double> y;
        vector<double> z;
        vector<int> xint;
        vector<int> yint;
        vector<int> zint;
	vector<double> E;
	int xpix,ypix,zpix;
	double xdist,ydist,zdist;
	cout<<"Start "<<centr_dist_x2<<"\t"<<centr_dist_y2<<"\t"<<centr_dist_z2<<endl;
        for(int xi=0;xi<140;xi++) {
                for(int yi=0;yi<128;yi++) {
                        for(int zi=0;zi<timetemp;zi++) {
//                              if((hcurrenttrack3D->GetBinContent(x,y,z)>0)) {//choose hbeamtrack for beam and hcurrenttrack3D for decay
                                if((hcurrenttrack3D->GetBinContent(xi+1,yi+1,zi+1)>0)) {//choose hbeamtrack for beam and hcurrenttrack3D for decay
                                        double maxE;
                                        maxE=hcurrenttrack3D->GetBinContent(xi+1,yi+1,zi+1);
                                        xpix=xi;
                                        ypix=yi;
                                        zpix=zi;
					for(int i=0;i<10;i++) {
	                                        Dist(xpix,ypix,zpix,xdist,ydist,zdist);
						x.push_back(xdist);
						y.push_back(ydist);
						z.push_back(zdist);
						xint.push_back(xi);
						yint.push_back(yi);
						zint.push_back(zi);
						E.push_back(maxE);
						if(good_evt<draw_events) {
							hdistxy1[good_evt]->Fill(xdist,ydist);
							hdistyz1[good_evt]->Fill(ydist,zdist);
							hdistxz1[good_evt]->Fill(xdist,zdist);
//							hdistxy2[good_evt]->Fill(xdist,ydist);
//							hdistyz2[good_evt]->Fill(ydist,zdist);
//							hdistxz2[good_evt]->Fill(xdist,zdist);
						}
					}
				}
			}
		}
	}
//
        const int npoints=x.size();
        double x_[npoints],y_[npoints],z_[npoints];
        for(int i=0;i<npoints;i++) {
                x_[i]=x.at(i);
                y_[i]=y.at(i);
                z_[i]=z.at(i);
        }
	if(npoints<3) {
		if(verbose) cout<<"Cannot RANSAC with "<<npoints<<" points"<<endl;
		return;
	}
//      RANSAC
        const int iterations = 100;
        const double threshold=3;
        int bestcounter=0;
        int bestselect1=0,bestselect2=0;
        for(int i=0;i<iterations;i++) {
                vector<int> pointid;
                int counter=0;
                int select1=int(npoints*rnd3->Rndm());
//                int select2=int(npoints*rnd3->Rndm());
//                double V1[3]={x.at(select1),y.at(select1),z.at(select1)},V2[3]={x.at(select2),y.at(select2),z.at(select2)};//points
                double V1[3]={x.at(select1),y.at(select1),z.at(select1)},V2[3]={centr_dist_x2,centr_dist_y2,centr_dist_z2};//points
                double L[3]={x.at(select1)-centr_dist_x2,y.at(select1)-centr_dist_y2,z.at(select1)-centr_dist_z2};
                double linelength=sqrt(L[0]*L[0]+L[1]*L[1]+L[2]*L[2]);
                for(int m=0;m<3;m++) L[m]/=linelength;//normalize the connecting line
//      Distance between point P and line r(t)=V1+t*L = d(P,L)=|(P-V1) x L|/|L|;
                for(int j=0;j<npoints;j++) {
			double best_dist=0;
                        double P[3]={x.at(j),y.at(j),z.at(j)};
                        double diffV[3]={0.};
                        for(int m=0;m<3;m++) diffV[m]=P[m]-V1[m];
//              diffV x L = i(diff_y*L_z-diff_z*L_y) - j(diff_x*L_z-diff_z*L_y) + k(diff_x*L_y-diff_y*L_x)
                        double cross[3]={0.};
                        cross[0]=diffV[1]*L[2]-diffV[2]*L[1];
                        cross[1]=-diffV[0]*L[2]+diffV[2]*L[0];
                        cross[2]=diffV[0]*L[1]-diffV[1]*L[0];
                        double distance=sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
//			cout<<P[0]<<"\t"<<P[1]<<"\t"<<P[2]<<"\t"<<distance<<"\t"<<V1[0]<<"\t"<<V1[1]<<"\t"<<V1[2]<<"\t"<<V2[0]<<"\t"<<V2[1]<<"\t"<<V2[2]<<endl;

                        if(distance<threshold) {
                                counter++;
                                pointid.push_back(j);
                        }
                }
                if(counter>bestcounter) {
//                        if(verbose)cout<<counter<<" out of "<<npoints<<endl;
                        bestcounter=counter;
                        bestselect1=select1;
//                        bestselect2=select2;
                }
        }
	

        double xbest[2],ybest[2],zbest[2];
        xbest[0]=x.at(bestselect1);
        xbest[1]=centr_dist_x2;
        ybest[0]=y.at(bestselect1);
        ybest[1]=centr_dist_y2;
        zbest[0]=z.at(bestselect1);
        zbest[1]=centr_dist_z2;
        double t[3]={xbest[0]-xbest[1],ybest[0]-ybest[1],zbest[0]-zbest[1]};
        double linex[1000],liney[1000],linez[1000];
        for(int i=0;i<50;i++) {
                linex[i]=xbest[0]+2*(i)*t[0];
                liney[i]=ybest[0]+2*(i)*t[1];
                linez[i]=zbest[0]+2*(i)*t[2];
		if(good_evt<draw_events) {
			hdistxy1[good_evt]->Fill(linex[i],liney[i],1e3);
			hdistyz1[good_evt]->Fill(liney[i],linez[i],1e3);
			hdistxz1[good_evt]->Fill(linex[i],linez[i],1e3);
		}
        }
	double totaldx=xbest[0]-xbest[1];
	double totaldy=ybest[0]-ybest[1];
	double totaldz=zbest[0]-zbest[1];
	double maxlength=sqrt((xbest[0]-xbest[1])*(xbest[0]-xbest[1])+(ybest[0]-ybest[1])*(ybest[0]-ybest[1])+(zbest[0]-zbest[1])*(zbest[0]-zbest[1]));
        length[0]=maxlength;
        palpha[0][0]=totaldx;
        palpha[0][1]=totaldy;
        palpha[0][2]=totaldz;
	scatterangle=atan(sqrt(totaldx*totaldx+totaldz*totaldz)/totaldy);
	chi_linear=0.;
	double total_w_energy=0.;
	double best_dist=0;

	for(int j=0;j<npoints;j++) {
	                double V1[3]={x.at(bestselect1),y.at(bestselect1),z.at(bestselect1)},V2[3]={centr_dist_x2,centr_dist_y2,centr_dist_z2};
	                double L[3]={x.at(bestselect1)-centr_dist_x2,y.at(bestselect1)-centr_dist_y2,z.at(bestselect1)-centr_dist_z2};
	                double linelength=sqrt(L[0]*L[0]+L[1]*L[1]+L[2]*L[2]);
        	        for(int m=0;m<3;m++) L[m]/=linelength;//normalize the connecting line
//			cout<<sqrt(L[0]*L[0]+L[1]*L[1]+L[2]*L[2])<<endl;
                        double P[3]={x.at(j),y.at(j),z.at(j)};
                        double diffV[3]={0.};
                        for(int m=0;m<3;m++) diffV[m]=P[m]-V2[m];
//              diffV x L = i(diff_y*L_z-diff_z*L_y) - j(diff_x*L_z-diff_z*L_y) + k(diff_x*L_y-diff_y*L_x)
                        double cross[3]={0.};
                        cross[0]=diffV[1]*L[2]-diffV[2]*L[1];
                        cross[1]=-diffV[0]*L[2]+diffV[2]*L[0];
                        cross[2]=diffV[0]*L[1]-diffV[1]*L[0];
                        double distance=sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
			double distance_from_vertex=0.;
			distance_from_vertex+=TMath::Power(centr_dist_x2-P[0],2);
			distance_from_vertex+=TMath::Power(centr_dist_y2-P[1],2);
			distance_from_vertex+=TMath::Power(centr_dist_z2-P[2],2);
			distance_from_vertex=sqrt(distance_from_vertex);
//			cout<<P[0]<<"\t"<<P[1]<<"\t"<<P[2]<<"\t"<<distance_from_vertex<<"\t"<<V1[0]<<"\t"<<V1[1]<<"\t"<<V1[2]<<"\t"<<V2[0]<<"\t"<<V2[1]<<"\t"<<V2[2]<<"\t"<<centr_dist_x2<<"\t"<<centr_dist_y2<<"\t"<<centr_dist_z2<<endl;
//			cout<<xint.at(j)<<"\t"<<yint.at(j)<<"\t"<<zint.at(j)<<endl;
//			if(distance<threshold) hcurrenttrack3D->SetBinContent(xint.at(j)+1,yint.at(j)+1,zint.at(j)+1,0);
			chi_linear+=distance*distance*E.at(j);
			total_w_energy+=E.at(j);
			if(distance_from_vertex>best_dist) {
				for(int m=0;m<3;m++) arm_point[m]=P[m];
//				cout<<"new best "<<best_dist<<endl;

				best_dist=distance_from_vertex;
			}
	}
	chi_linear/=total_w_energy;
	if(good_evt<draw_events) {
		hdistxy1[good_evt]->Fill(arm_point[0],arm_point[1],1e5);
		hdistyz1[good_evt]->Fill(arm_point[1],arm_point[2],1e5);
		hdistxz1[good_evt]->Fill(arm_point[0],arm_point[2],1e5);
		hdistxy1[good_evt]->Fill(centr_dist_x2,centr_dist_y2,2e5);
		hdistyz1[good_evt]->Fill(centr_dist_y2,centr_dist_z2,2e5);
		hdistxz1[good_evt]->Fill(centr_dist_x2,centr_dist_z2,2e5);
//		cout<<"end "<<arm_point[0]<<"\t"<<arm_point[1]<<"\t"<<arm_point[2]<<endl;
//		cout<<"start "<<centr_dist_x2<<"\t"<<centr_dist_y2<<"\t"<<centr_dist_z2<<endl;
	}

	return;
}
void Analysis::Ransac2() {
//	Get poivnt cloud
        vector<double> x;
        vector<double> y;
        vector<double> z;
        vector<int> xint;
        vector<int> yint;
        vector<int> zint;
	vector<double> E;
	int xpix,ypix,zpix;
	double xdist,ydist,zdist;
//	for(int i=0;i<5;i++) hcurrenttrack3D->Fill(rnd3->Rndm()*140,rnd3->Rndm()*128,rnd3->Rndm()*512.);
        for(int xi=0;xi<140;xi++) {
                for(int yi=0;yi<128;yi++) {
                        for(int zi=0;zi<timetemp;zi++) {
//                              if((hcurrenttrack3D->GetBinContent(x,y,z)>0)) {//choose hbeamtrack for beam and hcurrenttrack3D for decay
                                if((hcurrenttrack3D->GetBinContent(xi+1,yi+1,zi+1)>0)) {//choose hbeamtrack for beam and hcurrenttrack3D for decay
                                        double maxE;
                                        maxE=hcurrenttrack3D->GetBinContent(xi+1,yi+1,zi+1);
                                        xpix=xi;
                                        ypix=yi;
                                        zpix=zi;
					for(int k=0;k<nsamples;k++) {
	                                        Dist(xpix,ypix,zpix,xdist,ydist,zdist);
						x.push_back(xdist);
						y.push_back(ydist);
						z.push_back(zdist);
	                                        pointx.push_back(xdist);
	                                        pointy.push_back(ydist);
	                                        pointz.push_back(zdist);
	                                        pointE.push_back(maxE);
						E.push_back(maxE);
						xint.push_back(xi);
						yint.push_back(yi);
						zint.push_back(zi);
						if(good_evt<draw_events) {
							hdistxy[good_evt]->Fill(xdist,ydist);
							hdistyz[good_evt]->Fill(ydist,zdist);
							hdistxz[good_evt]->Fill(xdist,zdist);
							hdistxy2[good_evt]->Fill(xdist,ydist);
							hdistyz2[good_evt]->Fill(ydist,zdist);
							hdistxz2[good_evt]->Fill(xdist,zdist);
							hdistxy2a[good_evt]->Fill(xdist,ydist);
							hdistyz2a[good_evt]->Fill(ydist,zdist);
							hdistxz2a[good_evt]->Fill(xdist,zdist);
							hdistxy2b[good_evt]->Fill(xdist,ydist);
							hdistyz2b[good_evt]->Fill(ydist,zdist);
							hdistxz2b[good_evt]->Fill(xdist,zdist);
							hdistxy2c[good_evt]->Fill(xdist,ydist);
							hdistyz2c[good_evt]->Fill(ydist,zdist);
							hdistxz2c[good_evt]->Fill(xdist,zdist);
						}
					}
				}
			}
		}
	}
//
        const int npoints=x.size();
        double x_[npoints],y_[npoints],z_[npoints];
        for(int i=0;i<npoints;i++) {
                x_[i]=x.at(i);
                y_[i]=y.at(i);
                z_[i]=z.at(i);
        }
	if(npoints<4) return;
//      RANSAC
        const int iterations = 200*npoints*npoints;
        const double threshold=10;
	double bestdistance=1e40;
        int bestselect1=0,bestselect2=0,bestselect3=0;
	double total_w_energy=0.;
        for(int i=0;i<iterations;i++) {
                vector<int> pointid;
		total_w_energy=0.;
                int counter=0;
                int select1=int(npoints*rnd3->Rndm());
                int select2=int(npoints*rnd3->Rndm());
                int select3=int(npoints*rnd3->Rndm());
		if(select1==select2 || select1==select3 || select2==select3) continue;//stop double selection
                double V1[3]={x.at(select1),y.at(select1),z.at(select1)};//vertex
//		Check the vertex is at least VAGUELY central to correspond to the correct neutron interaction location 
		if(fabs(V1[0]-122.5)>10.) {
//			i--;
			continue;
		}
		double V2[3]={x.at(select2),y.at(select2),z.at(select2)};//points
		double V3[3]={x.at(select3),y.at(select3),z.at(select3)};//points
//		Check to make sure we have two tracks going either side of the vertex in x-direction and V2 is the left one
		if(V2[0]>V1[0] || V3[0]<V1[0]) {
//			i--;
			continue;//do not have one track either side of vertex so retry a sample
		}
		if((V2[2]>V1[2] && V3[2]>V1[2]) || (V2[2]<V1[2] && V3[2]<V1[2])) {
//			i--;
			continue;//do not have one track either side of vertex so retry a sample
		}
                double L1[3]={x.at(select1)-x.at(select2),y.at(select1)-y.at(select2),z.at(select1)-z.at(select2)};
                double linelength1=sqrt(L1[0]*L1[0]+L1[1]*L1[1]+L1[2]*L1[2]);
                double L2[3]={x.at(select1)-x.at(select3),y.at(select1)-y.at(select3),z.at(select1)-z.at(select3)};
                double linelength2=sqrt(L2[0]*L2[0]+L2[1]*L2[1]+L2[2]*L2[2]);
                for(int m=0;m<3;m++) L1[m]/=linelength1;//normalize the connecting line
                for(int m=0;m<3;m++) L2[m]/=linelength2;//normalize the connecting line
//      Distance between point P and line r(t)=V1+t*L = d(P,L)=|(P-V1) x L|/|L|;
		double totaldistance=0;
                for(int j=0;j<npoints;j++) {
			double best_dist=0;
                        double P[3]={x.at(j),y.at(j),z.at(j)};
                        double diffV1[3]={0.},diffV2[3]={0.};
                        for(int m=0;m<3;m++) {
				diffV1[m]=P[m]-V2[m];
				diffV2[m]=P[m]-V3[m];
			}
//              diffV x L = i(diff_y*L_z-diff_z*L_y) - j(diff_x*L_z-diff_z*L_y) + k(diff_x*L_y-diff_y*L_x)
                        double cross1[3]={0.};
                        cross1[0]=diffV1[1]*L1[2]-diffV1[2]*L1[1];
                        cross1[1]=-diffV1[0]*L1[2]+diffV1[2]*L1[0];
                        cross1[2]=diffV1[0]*L1[1]-diffV1[1]*L1[0];
                        double cross2[3]={0.};
                        cross2[0]=diffV2[1]*L2[2]-diffV2[2]*L2[1];
                        cross2[1]=-diffV2[0]*L2[2]+diffV2[2]*L2[0];
                        cross2[2]=diffV2[0]*L2[1]-diffV2[1]*L2[0];
                        double distance1=sqrt(cross1[0]*cross1[0]+cross1[1]*cross1[1]+cross1[2]*cross1[2]);
                        double distance2=sqrt(cross2[0]*cross2[0]+cross2[1]*cross2[1]+cross2[2]*cross2[2]);
			if(distance1<distance2) totaldistance+=distance1*distance1*E.at(j);
			if(distance2<distance1) totaldistance+=distance2*distance2*E.at(j);
			total_w_energy+=E.at(j);
                }
		totaldistance/=(linelength1+linelength2);
                if(totaldistance<bestdistance) {
                        bestdistance=totaldistance;
                        bestselect1=select1;//vertex
                        bestselect2=select2;
			bestselect3=select3;
//			cout<<"BEST n,a: "<<x.at(bestselect1)<<"\t"<<y.at(bestselect1)<<"\t"<<z.at(bestselect1)<<endl;

                }
        }
	for(int k=0;k<2;k++) {//loop arms
		int thisselect=0;
		if(k==0) thisselect=bestselect2;
		if(k==1) thisselect=bestselect3;
	        double xbest[2],ybest[2],zbest[2];
	       	xbest[0]=x.at(bestselect1);
	        xbest[1]=x.at(thisselect);
	        ybest[0]=y.at(bestselect1);
	        ybest[1]=y.at(thisselect);
	        zbest[0]=z.at(bestselect1);
	        zbest[1]=z.at(thisselect);
	        double t[3]={xbest[0]-xbest[1],ybest[0]-ybest[1],zbest[0]-zbest[1]};
	        double linex[1000],liney[1000],linez[1000];
	        for(int i=0;i<1000;i++) {
	                linex[i]=xbest[0]+.1*(i-50.)*t[0];
	                liney[i]=ybest[0]+.1*(i-50.)*t[1];
	                linez[i]=zbest[0]+.1*(i-50.)*t[2];
			if(good_evt<draw_events) {
//				hdistxy2[good_evt]->Fill(linex[i],liney[i],1e3);
//				hdistyz2[good_evt]->Fill(liney[i],linez[i],1e3);
//				hdistxz2[good_evt]->Fill(linex[i],linez[i],1e3);
			}
		}
        }
	double best_dist1=0,best_dist2=0.,best_dist3=0.;
	for(int k=0;k<3;k++) arm_energy[k]=0.;
	for(int j=0;j<npoints;j++) {
	                double V1[3]={x.at(bestselect1),y.at(bestselect1),z.at(bestselect1)};//vertex
	                double V2[3]={x.at(bestselect2),y.at(bestselect2),z.at(bestselect2)};
	                double V3[3]={x.at(bestselect3),y.at(bestselect3),z.at(bestselect3)};
	                double L1[3]={x.at(bestselect1)-x.at(bestselect2),y.at(bestselect1)-y.at(bestselect2),z.at(bestselect1)-z.at(bestselect2)};
	                double linelength1=sqrt(L1[0]*L1[0]+L1[1]*L1[1]+L1[2]*L1[2]);
	                double L2[3]={x.at(bestselect1)-x.at(bestselect3),y.at(bestselect1)-y.at(bestselect3),z.at(bestselect1)-z.at(bestselect3)};
	                double linelength2=sqrt(L2[0]*L2[0]+L2[1]*L2[1]+L2[2]*L2[2]);
	                for(int m=0;m<3;m++) L1[m]/=linelength1;//normalize the connecting line
	                for(int m=0;m<3;m++) L2[m]/=linelength2;//normalize the connecting line
                        double P[3]={x.at(j),y.at(j),z.at(j)};
                        double diffV1[3]={0.},diffV2[3]={0.};
                        for(int m=0;m<3;m++) {
				diffV1[m]=P[m]-V2[m];
				diffV2[m]=P[m]-V3[m];
			}
//              diffV x L = i(diff_y*L_z-diff_z*L_y) - j(diff_x*L_z-diff_z*L_y) + k(diff_x*L_y-diff_y*L_x)
                        double cross1[3]={0.};
                        cross1[0]=diffV1[1]*L1[2]-diffV1[2]*L1[1];
                        cross1[1]=-diffV1[0]*L1[2]+diffV1[2]*L1[0];
                        cross1[2]=diffV1[0]*L1[1]-diffV1[1]*L1[0];
                        double cross2[3]={0.};
                        cross2[0]=diffV2[1]*L2[2]-diffV2[2]*L2[1];
                        cross2[1]=-diffV2[0]*L2[2]+diffV2[2]*L2[0];
                        cross2[2]=diffV2[0]*L2[1]-diffV2[1]*L2[0];
                        double distance1=sqrt(cross1[0]*cross1[0]+cross1[1]*cross1[1]+cross1[2]*cross1[2]);
                        double distance2=sqrt(cross2[0]*cross2[0]+cross2[1]*cross2[1]+cross2[2]*cross2[2]);
			if(distance1<distance2) {
				arm_energy[0]+=E.at(j);
			}
			else {
				arm_energy[1]+=E.at(j);
			}
			double distance_from_vertex=0.;
			distance_from_vertex+=TMath::Power(V1[0]-P[0],2);
			distance_from_vertex+=TMath::Power(V1[1]-P[1],2);
			distance_from_vertex+=TMath::Power(V1[2]-P[2],2);
			distance_from_vertex=sqrt(distance_from_vertex);
//			cout<<P[0]<<"\t"<<P[1]<<"\t"<<P[2]<<"\t\t"<<distance_from_vertex<<"\t"<<distance1<<"\t"<<best_dist1<<"\t"<<distance2<<"\t"<<best_dist2<<endl;
			if(distance_from_vertex>best_dist1 && distance1<threshold && P[0]<V1[0]) {//New furthest away and on the line
				for(int m=0;m<3;m++) arm_point1[m]=P[m];
				arm_pix1[0]=xint.at(j);
				arm_pix1[1]=yint.at(j);
				arm_pix1[2]=zint.at(j);
				best_dist1=distance_from_vertex;
//				cout<<best_dist1<<"\t"<<distance_from_vertex<<endl;
				length[0]=best_dist1;
				palpha[0][0]=centr_dist_x-P[0];
				palpha[0][0]=centr_dist_y-P[1];
				palpha[0][0]=centr_dist_z-P[2];
				double totaldx=-(centr_dist_x-P[0]);
				double totaldy=-(centr_dist_y-P[1]);
				double totaldz=-(centr_dist_z-P[2]);
				scattering_theta[0]=atan2(sqrt(totaldx*totaldx+totaldz*totaldz),totaldy)/d2r;
//				cout<<"Scattering 0 "<<totaldx<<"\t"<<totaldy<<"\t"<<totaldz<<"\t"<<scattering_theta[0]<<endl;
				for(int m=0;m<3;m++)endpoint[0][m]=P[m];
			}
			if(distance_from_vertex>best_dist2 && distance2<threshold && P[0]>V1[0]) {
				for(int m=0;m<3;m++) arm_point2[m]=P[m];
				best_dist2=distance_from_vertex;
				arm_pix2[0]=xint.at(j);
				arm_pix2[1]=yint.at(j);
				arm_pix2[2]=zint.at(j);

				length[1]=best_dist2;
				palpha[1][0]=centr_dist_x-P[0];
				palpha[1][1]=centr_dist_y-P[1];
				palpha[1][2]=centr_dist_z-P[2];
				double totaldx=-(centr_dist_x-P[0]);
				double totaldy=-(centr_dist_y-P[1]);
				double totaldz=-(centr_dist_z-P[2]);
				scattering_theta[1]=atan2(sqrt(totaldx*totaldx+totaldz*totaldz),totaldy)/d2r;
//				cout<<"Scattering 1 "<<totaldx<<"\t"<<totaldy<<"\t"<<totaldz<<"\t"<<scattering_theta[1]<<endl;
				for(int m=0;m<3;m++)endpoint[1][m]=P[m];
			}
	}
	centr_dist_x=x.at(bestselect1);
	centr_dist_y=y.at(bestselect1);
	centr_dist_z=z.at(bestselect1);
//	cout<<centr_dist_x<<"\t"<<endpoint[0][0]<<"\t"<<endpoint[1][0]<<endl;
	if(good_evt<draw_events) {
		hdistxy2[good_evt]->Fill(endpoint[0][0],endpoint[0][1],1e8);
		hdistxy2[good_evt]->Fill(endpoint[1][0],endpoint[1][1],1e8);
		hdistxz2[good_evt]->Fill(endpoint[0][0],endpoint[0][2],1e8);
		hdistxz2[good_evt]->Fill(endpoint[1][0],endpoint[1][2],1e8);
		hdistyz2[good_evt]->Fill(endpoint[0][1],endpoint[0][2],1e8);
		hdistyz2[good_evt]->Fill(endpoint[1][1],endpoint[1][2],1e8);
		hdistxy2[good_evt]->Fill(x.at(bestselect1),y.at(bestselect1),3e8);
		hdistyz2[good_evt]->Fill(y.at(bestselect1),z.at(bestselect1),3e8);
		hdistxz2[good_evt]->Fill(x.at(bestselect1),z.at(bestselect1),3e8);
//		hdistxy[good_evt]->Fill(endpoint[0][0],endpoint[0][1],1e8);
//		hdistxy[good_evt]->Fill(endpoint[1][0],endpoint[1][1],1e8);
//		hdistxz[good_evt]->Fill(endpoint[0][0],endpoint[0][2],1e8);
//		hdistxz[good_evt]->Fill(endpoint[1][0],endpoint[1][2],1e8);
//		hdistyz[good_evt]->Fill(endpoint[0][1],endpoint[0][2],1e8);
//		hdistyz[good_evt]->Fill(endpoint[1][1],endpoint[1][2],1e8);
//		hdistxy[good_evt]->Fill(x.at(bestselect1),y.at(bestselect1),3e8);
//		hdistyz[good_evt]->Fill(y.at(bestselect1),z.at(bestselect1),3e8);
//		hdistxz[good_evt]->Fill(x.at(bestselect1),z.at(bestselect1),3e8);
//		cout<<x.at(bestselect1)<<"\t"<<y.at(bestselect1)<<"\t"<<z.at(bestselect1)<<endl;
//		htrackxy[good_evt]->Fill(arm_pix1[0],arm_pix1[1],1e6);
//		htrackyz[good_evt]->Fill(arm_pix1[1],arm_pix1[2],1e6);
//		htrackxz[good_evt]->Fill(arm_pix1[0],arm_pix1[2],1e6);
//		htrackxy[good_evt]->Fill(arm_pix2[0],arm_pix2[1],1e6);
//		htrackyz[good_evt]->Fill(arm_pix2[1],arm_pix2[2],1e6);
//		htrackxz[good_evt]->Fill(arm_pix2[0],arm_pix2[2],1e6);
//		htrackxy[good_evt]->Fill(xint.at(bestselect1),yint.at(bestselect1),1e6);
//		htrackyz[good_evt]->Fill(yint.at(bestselect1),zint.at(bestselect1),1e6);
//		htrackxz[good_evt]->Fill(xint.at(bestselect1),zint.at(bestselect1),1e6);
	}

	return;
}
void Analysis::Ransac3() {
//	Get point cloud
	chisqrdHoyle=0.;
        vector<double> x;
        vector<double> y;
        vector<double> z;
        vector<int> xint;
        vector<int> yint;
        vector<int> zint;
	vector<double> E;
	int xpix,ypix,zpix;
	int selend;
	double maxy;
	double xdist,ydist,zdist;
//	for(int i=0;i<5;i++) hcurrenttrack3D->Fill(rnd3->Rndm()*140,rnd3->Rndm()*128,rnd3->Rndm()*512.);
        for(int xi=0;xi<140;xi++) {
                for(int yi=0;yi<128;yi++) {
                        for(int zi=0;zi<timetemp;zi++) {
//                              if((hcurrenttrack3D->GetBinContent(x,y,z)>0)) {//choose hbeamtrack for beam and hcurrenttrack3D for decay
                                if((hcurrenttrack3D->GetBinContent(xi+1,yi+1,zi+1)>0)) {//choose hbeamtrack for beam and hcurrenttrack3D for decay
                                        double maxE;
                                        maxE=hcurrenttrack3D->GetBinContent(xi+1,yi+1,zi+1);
                                        xpix=xi;
                                        ypix=yi;
                                        zpix=zi;
                                        Dist(xpix,ypix,zpix,xdist,ydist,zdist);
					x.push_back(xdist);
					y.push_back(ydist);
					z.push_back(zdist);
					E.push_back(maxE);
					xint.push_back(xi);
					yint.push_back(yi);
					zint.push_back(zi);
					if(ydist>maxy) {
						maxy=ydist;
						selend=x.size()-1;
					}
					if(good_evt<draw_events) {
						hdistxy[good_evt]->Fill(xdist,ydist);
						hdistyz[good_evt]->Fill(ydist,zdist);
						hdistxz[good_evt]->Fill(xdist,zdist);
						hdistxy3[good_evt]->Fill(xdist,ydist);
						hdistyz3[good_evt]->Fill(ydist,zdist);
						hdistxz3[good_evt]->Fill(xdist,zdist);
//						hdistxy1[good_evt]->Fill(xdist,ydist);
//						hdistyz1[good_evt]->Fill(ydist,zdist);
//						hdistxz1[good_evt]->Fill(xdist,zdist);

					}
				}
			}
		}
	}
//
        const int npoints=x.size();
        double x_[npoints],y_[npoints],z_[npoints];
        for(int i=0;i<npoints;i++) {
                x_[i]=x.at(i);
                y_[i]=y.at(i);
                z_[i]=z.at(i);
        }
	if(npoints<4) return;
//      RANSAC
//        const int iterations = 10000;
        const int iterations = npoints*npoints;
        const double threshold=5;
	double bestdistance=1e40;
        int bestselect1=0,bestselect2=0,bestselect3=0;
	double total_w_energy=0.;
	double rcentr_dist_x=centr_dist_x;
	double rcentr_dist_y=centr_dist_y;
	double rcentr_dist_z=centr_dist_z;
        for(int i=0;i<iterations;i++) {
                vector<int> pointid;
		total_w_energy=0.;
                int counter=0;
                int select1=int(npoints*rnd3->Rndm());
                int select2=int(npoints*rnd3->Rndm());
                int select3=int(npoints*rnd3->Rndm());
//		rcentr_dist_x=x.at(select1);
//		rcentr_dist_y=y.at(select1);
//		rcentr_dist_z=z.at(select1);
		select1=selend;
                select2=floor(i/double(npoints));
                select3=i%npoints;

		if(select1==select2 || select1==select3 ||select2==select3) continue;
                double V1[3]={x.at(select1),y.at(select1),z.at(select1)};
		double V2[3]={x.at(select2),y.at(select2),z.at(select2)};//points
		double V3[3]={x.at(select3),y.at(select3),z.at(select3)};//points
                double L1[3]={x.at(select1)-rcentr_dist_x,y.at(select1)-rcentr_dist_y,z.at(select1)-rcentr_dist_z};
                double linelength1=sqrt(L1[0]*L1[0]+L1[1]*L1[1]+L1[2]*L1[2]);
                double L2[3]={x.at(select2)-rcentr_dist_x,y.at(select2)-rcentr_dist_y,z.at(select2)-rcentr_dist_z};
                double linelength2=sqrt(L2[0]*L2[0]+L2[1]*L2[1]+L2[2]*L2[2]);
                double L3[3]={x.at(select3)-rcentr_dist_x,y.at(select3)-rcentr_dist_y,z.at(select3)-rcentr_dist_z};
                double linelength3=sqrt(L3[0]*L3[0]+L3[1]*L3[1]+L3[2]*L3[2]);
                for(int m=0;m<3;m++) L1[m]/=linelength1;//normalize the connecting line
                for(int m=0;m<3;m++) L2[m]/=linelength2;//normalize the connecting line
                for(int m=0;m<3;m++) L3[m]/=linelength3;//normalize the connecting line
//      Distance between point P and line r(t)=V1+t*L = d(P,L)=|(P-V1) x L|/|L|;
		double totaldistance=0;
                for(int j=0;j<npoints;j++) {
			double best_dist=0;
                        double P[3]={x.at(j),y.at(j),z.at(j)};
                        double diffV1[3]={0.},diffV2[3]={0.},diffV3[3]={0.};
                        for(int m=0;m<3;m++) {
				diffV1[m]=P[m]-V1[m];
				diffV2[m]=P[m]-V2[m];
				diffV3[m]=P[m]-V3[m];
			}
//              diffV x L = i(diff_y*L_z-diff_z*L_y) - j(diff_x*L_z-diff_z*L_y) + k(diff_x*L_y-diff_y*L_x)
                        double cross1[3]={0.};
                        cross1[0]=diffV1[1]*L1[2]-diffV1[2]*L1[1];
                        cross1[1]=-diffV1[0]*L1[2]+diffV1[2]*L1[0];
                        cross1[2]=diffV1[0]*L1[1]-diffV1[1]*L1[0];
                        double cross2[3]={0.};
                        cross2[0]=diffV2[1]*L2[2]-diffV2[2]*L2[1];
                        cross2[1]=-diffV2[0]*L2[2]+diffV2[2]*L2[0];
                        cross2[2]=diffV2[0]*L2[1]-diffV2[1]*L2[0];
                        double cross3[3]={0.};
                        cross3[0]=diffV3[1]*L3[2]-diffV3[2]*L3[1];
                        cross3[1]=-diffV3[0]*L3[2]+diffV3[2]*L3[0];
                        cross3[2]=diffV3[0]*L3[1]-diffV3[1]*L3[0];
                        double distance1=sqrt(cross1[0]*cross1[0]+cross1[1]*cross1[1]+cross1[2]*cross1[2]);
                        double distance2=sqrt(cross2[0]*cross2[0]+cross2[1]*cross2[1]+cross2[2]*cross2[2]);
                        double distance3=sqrt(cross3[0]*cross3[0]+cross3[1]*cross3[1]+cross3[2]*cross3[2]);
//			cout<<"D: "<<distance1<<"\t"<<distance2<<"\t"<<distance3<<endl;
//			cout<<P[0]<<"\t"<<P[1]<<"\t"<<P[2]<<"\t"<<distance<<"\t"<<V1[0]<<"\t"<<V1[1]<<"\t"<<V1[2]<<"\t"<<V2[0]<<"\t"<<V2[1]<<"\t"<<V2[2]<<endl;
			if(distance1<distance2 && distance1<distance3) totaldistance+=distance1*distance1*E.at(j);
			if(distance2<distance1 && distance2<distance3) totaldistance+=distance2*distance2*E.at(j);
			if(distance3<distance2 && distance3<distance2) totaldistance+=distance3*distance3*E.at(j);
			total_w_energy+=E.at(j);
                }
                if(totaldistance<bestdistance) {
//			if(totaldistance>1e2) totaldistance=1e2;
                        bestdistance=totaldistance;
			chi_triple=bestdistance/(total_w_energy);
                        bestselect1=select1;
                        bestselect2=select2;
			bestselect3=select3;
//			cout<<select1<<"\t"<<select2<<"\t"<<select3<<"\t"<<bestdistance<<endl;
                }
        }
	
	for(int k=0;k<3;k++) {//loop arms
		int thisselect=0;
		if(k==0) thisselect=bestselect1;
		if(k==1) thisselect=bestselect2;
		if(k==2) thisselect=bestselect3;
	        double xbest[2],ybest[2],zbest[2];
	       	xbest[0]=x.at(thisselect);
	        xbest[1]=rcentr_dist_x;
	        ybest[0]=y.at(thisselect);
	        ybest[1]=rcentr_dist_y;
	        zbest[0]=z.at(thisselect);
	        zbest[1]=rcentr_dist_z;
	        double t[3]={xbest[0]-xbest[1],ybest[0]-ybest[1],zbest[0]-zbest[1]};
	        double linex[1000],liney[1000],linez[1000];
	        for(int i=0;i<1000;i++) {
	                linex[i]=xbest[0]+.1*(i-50.)*t[0];
	                liney[i]=ybest[0]+.1*(i-50.)*t[1];
	                linez[i]=zbest[0]+.1*(i-50.)*t[2];
			if(good_evt<draw_events && 0) {
				hdistxy3[good_evt]->Fill(linex[i],liney[i],1e3);
				hdistyz3[good_evt]->Fill(liney[i],linez[i],1e3);
				hdistxz3[good_evt]->Fill(linex[i],linez[i],1e3);
			}
		}
        }

	arm1_E=0.;
	arm3_E=0.;
	arm3_E=0.;
	double best_dist1=0,best_dist2=0.,best_dist3=0.;
	double arm_threshold=100;
	for(int j=0;j<npoints;j++) {
	                double V1[3]={x.at(bestselect1),y.at(bestselect1),z.at(bestselect1)};
	                double V2[3]={x.at(bestselect2),y.at(bestselect2),z.at(bestselect2)};
	                double V3[3]={x.at(bestselect3),y.at(bestselect3),z.at(bestselect3)};
	                double L1[3]={x.at(bestselect1)-rcentr_dist_x,y.at(bestselect1)-rcentr_dist_y,z.at(bestselect1)-rcentr_dist_z};
	                double linelength1=sqrt(L1[0]*L1[0]+L1[1]*L1[1]+L1[2]*L1[2]);
	                double L2[3]={x.at(bestselect2)-rcentr_dist_x,y.at(bestselect2)-rcentr_dist_y,z.at(bestselect2)-rcentr_dist_z};
	                double linelength2=sqrt(L2[0]*L2[0]+L2[1]*L2[1]+L2[2]*L2[2]);
	                double L3[3]={x.at(bestselect3)-rcentr_dist_x,y.at(bestselect3)-rcentr_dist_y,z.at(bestselect3)-rcentr_dist_z};
	                double linelength3=sqrt(L3[0]*L3[0]+L3[1]*L3[1]+L3[2]*L3[2]);
	                for(int m=0;m<3;m++) L1[m]/=linelength1;//normalize the connecting line
	                for(int m=0;m<3;m++) L2[m]/=linelength2;//normalize the connecting line
	                for(int m=0;m<3;m++) L3[m]/=linelength3;//normalize the connecting line
                        double P[3]={x.at(j),y.at(j),z.at(j)};
                        double diffV1[3]={0.},diffV2[3]={0.},diffV3[3]={0.};
                        for(int m=0;m<3;m++) {
				diffV1[m]=P[m]-V1[m];
				diffV2[m]=P[m]-V2[m];
				diffV3[m]=P[m]-V3[m];
			}
//              diffV x L = i(diff_y*L_z-diff_z*L_y) - j(diff_x*L_z-diff_z*L_y) + k(diff_x*L_y-diff_y*L_x)
                        double cross1[3]={0.};
                        cross1[0]=diffV1[1]*L1[2]-diffV1[2]*L1[1];
                        cross1[1]=-diffV1[0]*L1[2]+diffV1[2]*L1[0];
                        cross1[2]=diffV1[0]*L1[1]-diffV1[1]*L1[0];
                        double cross2[3]={0.};
                        cross2[0]=diffV2[1]*L2[2]-diffV2[2]*L2[1];
                        cross2[1]=-diffV2[0]*L2[2]+diffV2[2]*L2[0];
                        cross2[2]=diffV2[0]*L2[1]-diffV2[1]*L2[0];
                        double cross3[3]={0.};
                        cross3[0]=diffV3[1]*L3[2]-diffV3[2]*L3[1];
                        cross3[1]=-diffV3[0]*L3[2]+diffV3[2]*L3[0];
                        cross3[2]=diffV3[0]*L3[1]-diffV3[1]*L3[0];
                        double distance1=sqrt(cross1[0]*cross1[0]+cross1[1]*cross1[1]+cross1[2]*cross1[2]);
                        double distance2=sqrt(cross2[0]*cross2[0]+cross2[1]*cross2[1]+cross2[2]*cross2[2]);
                        double distance3=sqrt(cross3[0]*cross3[0]+cross3[1]*cross3[1]+cross3[2]*cross3[2]);
//
//			cout<<x.at(j)<<"\t"<<y.at(j)<<"\t"<<z.at(j)<<"\t"<<distance1<<"\t"<<distance2<<"\t"<<distance3<<endl;
			double distance_from_vertex=0.;
			distance_from_vertex+=TMath::Power(rcentr_dist_x-P[0],2);
			distance_from_vertex+=TMath::Power(rcentr_dist_y-P[1],2);
			distance_from_vertex+=TMath::Power(rcentr_dist_z-P[2],2);
			distance_from_vertex=sqrt(distance_from_vertex);

			if(distance1<distance2 && distance1<distance3 && distance1<arm_threshold) {
				arm1_E+=E.at(j);
				if(good_evt<draw_events) {
					harmxy[good_evt]->Fill(x.at(j),y.at(j),1);
					harmyz[good_evt]->Fill(y.at(j),z.at(j),1);
					harmxz[good_evt]->Fill(x.at(j),z.at(j),1);
				}
				if(distance_from_vertex>best_dist1) {
					for(int m=0;m<3;m++) arm_point1[m]=P[m];
					best_dist1=distance_from_vertex;
					length[0]=best_dist1;
					palpha[0][0]=rcentr_dist_x-P[0];
					palpha[0][1]=rcentr_dist_y-P[1];
					palpha[0][2]=rcentr_dist_z-P[2];
				}
			chisqrdHoyle+=distance1*distance1;
			}
			if(distance2<distance1 && distance2<distance3 && distance2<arm_threshold) {
				arm2_E+=E.at(j);
				if(good_evt<draw_events) {
					harmxy[good_evt]->Fill(x.at(j),y.at(j),10);
					harmyz[good_evt]->Fill(y.at(j),z.at(j),10);
					harmxz[good_evt]->Fill(x.at(j),z.at(j),10);
				}
				if(distance_from_vertex>best_dist2) {
					for(int m=0;m<3;m++) arm_point2[m]=P[m];
					best_dist2=distance_from_vertex;
					length[1]=best_dist2;
					palpha[1][0]=rcentr_dist_x-P[0];
					palpha[1][1]=rcentr_dist_y-P[1];
					palpha[1][2]=rcentr_dist_z-P[2];
				}
			chisqrdHoyle+=distance2*distance2;

			}
			if(distance3<distance2 && distance3<distance1 && distance3<arm_threshold) {
				arm3_E+=E.at(j);
				if(good_evt<draw_events) {
					harmxy[good_evt]->Fill(x.at(j),y.at(j),100);
					harmyz[good_evt]->Fill(y.at(j),z.at(j),100);
					harmxz[good_evt]->Fill(x.at(j),z.at(j),100);
				}
				if(distance_from_vertex>best_dist3) {
					for(int m=0;m<3;m++) arm_point3[m]=P[m];
					best_dist3=distance_from_vertex;
					length[2]=best_dist3;
					palpha[2][0]=rcentr_dist_x-P[0];
					palpha[2][1]=rcentr_dist_y-P[1];
					palpha[2][2]=rcentr_dist_z-P[2];
				}
			chisqrdHoyle+=distance3*distance3;
			}
			chisqrdHoyle/=double(npoints);
/*
			if(P[1]>centr_dist_y) {
				if(distance_from_vertex>best_dist1 && distance1<0.2*distance_from_vertex+threshold && distance1<distance2 && distance1<distance3) {//New furthest away and on the line
					for(int m=0;m<3;m++) arm_point1[m]=P[m];
					best_dist1=distance_from_vertex;
					length[0]=best_dist1;
					palpha[0][0]=centr_dist_x-P[0];
					palpha[0][0]=centr_dist_y-P[1];
					palpha[0][0]=centr_dist_z-P[2];
				}
				if(distance_from_vertex>best_dist2 && distance2<0.2*distance_from_vertex+threshold  && distance2<distance1 && distance2<distance3) {
					for(int m=0;m<3;m++) arm_point2[m]=P[m];
					best_dist2=distance_from_vertex;
					length[1]=best_dist2;
					palpha[1][0]=centr_dist_x-P[0];
					palpha[1][1]=centr_dist_y-P[1];
					palpha[1][2]=centr_dist_z-P[2];
				}
				if(distance_from_vertex>best_dist3 && distance3<0.2*distance_from_vertex+threshold && distance3<distance1 && distance3<distance2) {
					for(int m=0;m<3;m++) arm_point3[m]=P[m];
					best_dist3=distance_from_vertex;
					length[2]=best_dist3;
					palpha[2][0]=centr_dist_x-P[0];
					palpha[2][1]=centr_dist_y-P[1];
					palpha[2][2]=centr_dist_z-P[2];
				}
			}
*/
	}
	if(good_evt<draw_events) {
//		cout<<arm_point1[0]<<"\t"<<arm_point1[1]<<"\t"<<arm_point1[2]<<endl;
//		cout<<arm_point2[0]<<"\t"<<arm_point2[1]<<"\t"<<arm_point2[2]<<endl;
//		cout<<arm_point3[0]<<"\t"<<arm_point3[1]<<"\t"<<arm_point3[2]<<endl;
		hdistxy3[good_evt]->Fill(arm_point1[0],arm_point1[1],1e8);
		hdistyz3[good_evt]->Fill(arm_point1[1],arm_point1[2],1e8);
		hdistxz3[good_evt]->Fill(arm_point1[0],arm_point1[2],1e8);
//		cout<<arm_point1[0]<<"\t"<<arm_point1[1]<<"\t"<<arm_point1[2]<<endl;
		hdistxy3[good_evt]->Fill(arm_point2[0],arm_point2[1],1e8);
		hdistyz3[good_evt]->Fill(arm_point2[1],arm_point2[2],1e8);
		hdistxz3[good_evt]->Fill(arm_point2[0],arm_point2[2],1e8);
//		cout<<arm_point2[0]<<"\t"<<arm_point2[1]<<"\t"<<arm_point2[2]<<endl;
		hdistxy3[good_evt]->Fill(arm_point3[0],arm_point3[1],1e8);
		hdistyz3[good_evt]->Fill(arm_point3[1],arm_point3[2],1e8);
		hdistxz3[good_evt]->Fill(arm_point3[0],arm_point3[2],1e8);
//
		hdistxy3[good_evt]->Fill(rcentr_dist_x,rcentr_dist_y,2e8);
		hdistyz3[good_evt]->Fill(rcentr_dist_y,rcentr_dist_z,2e8);
		hdistxz3[good_evt]->Fill(rcentr_dist_x,rcentr_dist_z,2e8);
//
		hdistxy3[good_evt]->SetTitle(Form("pixel energies for arms: %g %g %g",arm1_E,arm2_E,arm3_E));
		hdistyz3[good_evt]->SetTitle(Form("Chi-sqrd Hoyle: %g",chisqrdHoyle));
//		cout<<arm_point3[0]<<"\t"<<arm_point3[1]<<"\t"<<arm_point3[2]<<endl;
	}

	return;
}



double Analysis::Fitnalpha() {
//Use TMinuit minimizer to fit the best arms
//Chi-sqrd minimization of distance to nearest arms for all voxels (weighted by amplitude)
        ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2","");//creater minuit2 minimizer
        ROOT::Math::Functor f(this,&Analysis::CalcChisqrdDist,9);//functor allows for this class function to be used instead of a standalone function as normally required
        min->SetFunction(f);//set this functor which points to the chisqrd calculator
//      Set initial parameters
        Double_t vstart[9] = {centr_dist_x,centr_dist_y,centr_dist_z,endpoint[0][0],endpoint[0][1],endpoint[0][2],endpoint[1][0],endpoint[1][1],endpoint[1][2]};
        Double_t vstep[9] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};
        if((verbose || debug)) {
                cout<<evtno<<" Pre:\t";
                for(int n=0;n<9;n++) {
                        cout<<vstart[n]<<"\t";
                }
                cout<<"\n"<<endl;
        }
//      EDM is expected distance from the minimum
//      Set names and starting values and steps sizes
	looper=0;
	min->SetVariableLimits(0,1.110,135);
	min->SetVariableStepSize(0,3.5);
	min->SetVariableStepSize(1,3.5);
	min->SetVariableStepSize(2,3.5);
	min->SetVariableStepSize(3,3.5);
	min->SetVariableStepSize(4,3.5);
	min->SetVariableStepSize(5,3.5);
	min->SetVariableStepSize(6,3.5);
	min->SetVariableStepSize(7,3.5);
	min->SetVariableStepSize(8,3.5);
	min->SetVariableStepSize(9,3.5);
        min->SetVariable(0,"Center x",vstart[0],vstep[0]);
        min->SetVariable(1,"Center y",vstart[1],vstep[1]);
        min->SetVariable(2,"Center z",vstart[2],vstep[2]);
        min->SetVariable(3,"Arm1 x",vstart[3],vstep[3]);
        min->SetVariable(4,"Arm1 y",vstart[4],vstep[4]);
        min->SetVariable(5,"Arm1 z",vstart[5],vstep[5]);
        min->SetVariable(6,"Arm2 x",vstart[6],vstep[6]);
        min->SetVariable(7,"Arm2 y",vstart[7],vstep[7]);
        min->SetVariable(8,"Arm2 z",vstart[8],vstep[8]);
        min->SetPrintLevel(verbose?1:0);
//      min->SetPrintLevel(1);
        min->Minimize();
        const double *param_loop0 = min->X();
        double chisqrd_loop0=min->MinValue()+chisqrd_nalpha;
//	chisqrd_loop0=chisqrd_nalpha;

//
	looper=1;
	min->SetVariableLimits(0,1.115,130);
        min->SetVariable(0,"Center x",vstart[0],vstep[0]);
        min->SetVariable(1,"Center y",vstart[1],vstep[1]);
        min->SetVariable(2,"Center z",vstart[2],vstep[2]);
        min->SetVariable(3,"Arm1 x",vstart[3],vstep[3]);
        min->SetVariable(4,"Arm1 y",vstart[4],vstep[4]);
        min->SetVariable(5,"Arm1 z",vstart[5],vstep[5]);
        min->SetVariable(6,"Arm2 x",vstart[6],vstep[6]);
        min->SetVariable(7,"Arm2 y",vstart[7],vstep[7]);
        min->SetVariable(8,"Arm2 z",vstart[8],vstep[8]);
        min->SetPrintLevel(verbose?1:0);
//      min->SetPrintLevel(1);
        min->Minimize();
        const double *param_loop1 = min->X();
        double chisqrd_loop1=min->MinValue()+chisqrd_nalpha;
//	chisqrd_loop1=chisqrd_nalpha;
	double chisqrd=1e6;
	double param[9];
	cout<<select_nalpha<<"\t"<<chisqrd_loop0<<"\t"<<chisqrd_loop1<<endl;
	if(chisqrd_loop0<chisqrd_loop1) {
		for(int l=0;l<9;l++)param[l]=param_loop0[l];
		chisqrd=chisqrd_loop0;
		fit_theta[select_nalpha-1][0]=ftheta[0][0];
		fit_theta[select_nalpha-1][1]=ftheta[0][1];
		dlength[select_nalpha-1][0]=dE_light_[0];
		dlength[select_nalpha-1][1]=dE_heavy_[0];
		for(int k=0;k<3;k++) p_vertex[k]=param_loop0[k];
		for(int k=0;k<3;k++) p_end1[k]=param_loop0[k+3];
		for(int k=0;k<3;k++) p_end2[k]=param_loop0[k+6];

	}
	if(chisqrd_loop1<chisqrd_loop0) {
		for(int l=0;l<9;l++)param[l]=param_loop1[l];
		chisqrd=chisqrd_loop1;
		fit_theta[select_nalpha-1][0]=ftheta[1][0];
		fit_theta[select_nalpha-1][1]=ftheta[1][1];
		dlength[select_nalpha-1][0]=dE_light_[1];
		dlength[select_nalpha-1][1]=dE_heavy_[1];
		for(int k=0;k<3;k++) p_vertex[k]=param_loop1[k];
		for(int k=0;k<3;k++) p_end1[k]=param_loop1[k+3];
		for(int k=0;k<3;k++) p_end2[k]=param_loop1[k+6];

	}
	p_chisqrd=chisqrd;

//
//	cout<<select_nalpha<<"\t"<<chisqrd<<endl;
//      double *covariancematrix = min->GetCovMatrix();//filled as cov[i*ndim + j]
//        min->Hesse();
//      cout<<"STATUS "<<min->CovMatrixStatus()<<"\t"<<chisqrd<<endl;
//        double covariancematrix[144];
  //      min->GetCovMatrix(covariancematrix);//filled as cov[i*ndim + j]
//	for(int i=0;i<9;i++) {
//		cout<<i<<"\t"<<vstart[i]<<"\t"<<param[i]<<endl;
//	}
	vertexy=param[1];
	if(good_evt<draw_events) {
		for(double l=0;l<100;l+=0.1) {
			if(select_nalpha==1) {
				hdistxy2a[good_evt]->Fill(param[0]+(l/100.)*(param[3]-param[0]),param[1]+(l/100.)*(param[4]-param[1]));
				hdistxy2a[good_evt]->Fill(param[0]+(l/100.)*(param[6]-param[0]),param[1]+(l/100.)*(param[7]-param[1]));
				hdistxy2a[good_evt]->Fill(param[0],param[1],1e2);
				hdistxy2a[good_evt]->Fill(param[3],param[4],1e3);
				hdistxy2a[good_evt]->Fill(param[6],param[7],1e3);
//
				hdistyz2a[good_evt]->Fill(param[1]+(l/100.)*(param[4]-param[1]),param[2]+(l/100.)*(param[5]-param[2]));
				hdistyz2a[good_evt]->Fill(param[1]+(l/100.)*(param[7]-param[1]),param[2]+(l/100.)*(param[8]-param[2]));
				hdistyz2a[good_evt]->Fill(param[1],param[2],1e2);
				hdistyz2a[good_evt]->Fill(param[4],param[5],1e3);
				hdistyz2a[good_evt]->Fill(param[7],param[8],1e3);
//
				hdistxz2a[good_evt]->Fill(param[0]+(l/100.)*(param[3]-param[0]),param[2]+(l/100.)*(param[5]-param[2]));
				hdistxz2a[good_evt]->Fill(param[0]+(l/100.)*(param[6]-param[0]),param[2]+(l/100.)*(param[8]-param[2]));
				hdistxz2a[good_evt]->Fill(param[0],param[2],1e2);
				hdistxz2a[good_evt]->Fill(param[3],param[5],1e3);
				hdistxz2a[good_evt]->Fill(param[6],param[8],1e3);
			}
			if(select_nalpha==2) {
				hdistxy2b[good_evt]->Fill(param[0]+(l/100.)*(param[3]-param[0]),param[1]+(l/100.)*(param[4]-param[1]));
				hdistxy2b[good_evt]->Fill(param[0]+(l/100.)*(param[6]-param[0]),param[1]+(l/100.)*(param[7]-param[1]));
				hdistxy2b[good_evt]->Fill(param[0],param[1],1e2);
				hdistxy2b[good_evt]->Fill(param[3],param[4],1e3);
				hdistxy2b[good_evt]->Fill(param[6],param[7],1e3);
//
				hdistyz2b[good_evt]->Fill(param[1]+(l/100.)*(param[4]-param[1]),param[2]+(l/100.)*(param[5]-param[2]));
				hdistyz2b[good_evt]->Fill(param[1]+(l/100.)*(param[7]-param[1]),param[2]+(l/100.)*(param[8]-param[2]));
				hdistyz2b[good_evt]->Fill(param[1],param[2],1e2);
				hdistyz2b[good_evt]->Fill(param[4],param[5],1e3);
				hdistyz2b[good_evt]->Fill(param[7],param[8],1e3);
//
				hdistxz2b[good_evt]->Fill(param[0]+(l/100.)*(param[3]-param[0]),param[2]+(l/100.)*(param[5]-param[2]));
				hdistxz2b[good_evt]->Fill(param[0]+(l/100.)*(param[6]-param[0]),param[2]+(l/100.)*(param[8]-param[2]));
				hdistxz2b[good_evt]->Fill(param[0],param[2],1e2);
				hdistxz2b[good_evt]->Fill(param[3],param[5],1e3);
				hdistxz2b[good_evt]->Fill(param[6],param[8],1e3);
			}
			if(select_nalpha==3) {
				hdistxy2c[good_evt]->Fill(param[0]+(l/100.)*(param[3]-param[0]),param[1]+(l/100.)*(param[4]-param[1]));
				hdistxy2c[good_evt]->Fill(param[0]+(l/100.)*(param[6]-param[0]),param[1]+(l/100.)*(param[7]-param[1]));
				hdistxy2c[good_evt]->Fill(param[0],param[1],1e2);
				hdistxy2c[good_evt]->Fill(param[3],param[4],1e3);
				hdistxy2c[good_evt]->Fill(param[6],param[7],1e3);
//
				hdistyz2c[good_evt]->Fill(param[1]+(l/100.)*(param[4]-param[1]),param[2]+(l/100.)*(param[5]-param[2]));
				hdistyz2c[good_evt]->Fill(param[1]+(l/100.)*(param[7]-param[1]),param[2]+(l/100.)*(param[8]-param[2]));
				hdistyz2c[good_evt]->Fill(param[1],param[2],1e2);
				hdistyz2c[good_evt]->Fill(param[4],param[5],1e3);
				hdistyz2c[good_evt]->Fill(param[7],param[8],1e3);
//
				hdistxz2c[good_evt]->Fill(param[0]+(l/100.)*(param[3]-param[0]),param[2]+(l/100.)*(param[5]-param[2]));
				hdistxz2c[good_evt]->Fill(param[0]+(l/100.)*(param[6]-param[0]),param[2]+(l/100.)*(param[8]-param[2]));
				hdistxz2c[good_evt]->Fill(param[0],param[2],1e2);
				hdistxz2c[good_evt]->Fill(param[3],param[5],1e3);
				hdistxz2c[good_evt]->Fill(param[6],param[8],1e3);
			}

		}
	}
//	chisqrd=chisqrd_nalpha;
	cout<<"Fit n alpha with chisqrd of "<<chisqrd<<endl;

	return chisqrd;
}


void Analysis::Cross(double x1,double y1,double z1,double x2,double y2,double z2,double &x3, double &y3, double &z3) {
//     Calculate cross product of (x1,y1,z1) with (x2,y2,z2)
        x3=y1*z2-z1*y2;
        y3=z1*x2-x1*z2;
        z3=x1*y2-x2*y1;
        return;
}

double Analysis::FitLine() {
       ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2","");//creater minuit2 minimizer
        ROOT::Math::Functor f(this,&Analysis::CalcChisqrdLine,6);//functor allows for this class function to be used instead of a standalone$
        min->SetFunction(f);//set this functor which points to the chisqrd calculator
//      Set initial parameters
        Double_t vstart[6] = {centr_dist_x,centr_dist_y,centr_dist_z,endpoint[0][0],endpoint[0][1],endpoint[0][2]};
        Double_t vstep[6] = {0.1,0.1,0.1,0.1,0.1,0.1};
        if((verbose || debug)) {
                cout<<evtno<<" Pre:\t";
                for(int n=0;n<6;n++) {
                        cout<<vstart[n]<<"\t";
                }
                cout<<"\n"<<endl;
        }
//      EDM is expected distance from the minimum
//      Set names and starting values and steps sizes
        min->SetVariable(0,"Center x",vstart[0],vstep[0]);
        min->SetVariable(1,"Center y",vstart[1],vstep[1]);
        min->SetVariable(2,"Center z",vstart[2],vstep[2]);
        min->SetVariable(3,"Arm1 x",vstart[3],vstep[3]);
        min->SetVariable(4,"Arm1 y",vstart[4],vstep[4]);
        min->SetVariable(5,"Arm1 z",vstart[5],vstep[5]);
        min->SetPrintLevel(verbose?1:0);
//      min->SetPrintLevel(1);
        min->Minimize();
        const double *param = min->X();
        const double chisqrd=min->MinValue();
	if(good_evt<draw_events) {
		for(double l=-5;l<5;l+=0.01) {
//			cout<<param[0]+l*(param[3]-param[0])<<"\t"<<param[1]+l*(param[4]-param[1])<<"\t"<<param[2]+l*(param[5]-param[2])<<endl;
			hdistxy1[good_evt]->Fill(param[0]+l*(param[3]-param[0]),param[1]+l*(param[4]-param[1]));
			hdistyz1[good_evt]->Fill(param[1]+l*(param[4]-param[1]),param[2]+l*(param[5]-param[2]));
			hdistxz1[good_evt]->Fill(param[0]+l*(param[3]-param[0]),param[2]+l*(param[5]-param[2]));
		}
	}
	double dx=fabs(param[3]-param[0]);
	double dy=fabs(param[4]-param[1]);
	double dz=fabs(param[5]-param[2]);
	linetheta=atan2(sqrt(dx*dx+dz*dz),dy)/d2r;
	linedist=sqrt(dx*dx+dy*dy+dz*dz);
	helscat->Fill(linetheta,linedist);
	return chisqrd;
};
double Analysis::CalcChisqrdLine(const double *p) {
	double chisqrd=0.;
        double L1[3]={p[3]-p[0],p[4]-p[1],p[5]-p[2]};//current line fit
        double linelength1=sqrt(L1[0]*L1[0]+L1[1]*L1[1]+L1[2]*L1[2]);
        for(int m=0;m<3;m++) L1[m]/=linelength1;//normalize the connecting line
	double minimum_y=1000;
	double maximum_y=0;
        for(unsigned int num=0;num<pointx.size();num++) {
                double cross1x,cross1y,cross1z,bestdiff=1000.;
//      Distance between point P and line r(t)=V1+t*L = d(P,L)=|(P-V1) x L|/|L|;
                double totaldistance=0;
                double P[3]={pointx.at(num),pointy.at(num),pointz.at(num)};//current point
		if(P[1]>maximum_y) maximum_y=P[1];
		if(P[1]<minimum_y) minimum_y=P[1];
                double diffV1[3]={0.};
                for(int m=0;m<3;m++) {
                        diffV1[m]=P[m]-p[m];
                }
//              diffV x L = i(diff_y*L_z-diff_z*L_y) - j(diff_x*L_z-diff_z*L_y) + k(diff_x*L_y-diff_y*L_x)
                double cross1[3]={0.};
                cross1[0]=diffV1[1]*L1[2]-diffV1[2]*L1[1];
                cross1[1]=-diffV1[0]*L1[2]+diffV1[2]*L1[0];
                cross1[2]=diffV1[0]*L1[1]-diffV1[1]*L1[0];
                double distance1=sqrt(cross1[0]*cross1[0]+cross1[1]*cross1[1]+cross1[2]*cross1[2]);
		chisqrd+=distance1*distance1;
	}
	chisqrd+=TMath::Power(p[1]-minimum_y,2.);//align to end of track
	chisqrd+=TMath::Power(p[4]-maximum_y,2.);//align to end of track
//	chisqrd-=linelength1*linelength1;
	double ndof=double(pointx.size())-6;
	chisqrd/=ndof;
	return chisqrd;
}

double Analysis::GetLabAngle12C(double cmangle) {
        cmangle*=d2r;
        double ma=1;
        double mA=12.;
        double mb=4.;
        double mB=9;
        double Q=-5.70205;
        double Ea=Ebeam;
        double gamma=sqrt(((ma*mb)/(mA*mB)) * ((Ea)/((1.+ma/mA)*Q+Ea)));
        double coscm=cos(cmangle);
        double coslab=(gamma+coscm)/sqrt(1+gamma*gamma+2*gamma*coscm);
        return acos(coslab)/d2r;

}

double Analysis::GetLabAngle16O(double cmangle) {
        cmangle*=d2r;
        double ma=1;
        double mA=16.;
        double mb=4.;
        double mB=13;
        double Q=-2.2156;
        double Ea=Ebeam;
        double gamma=sqrt(((ma*mb)/(mA*mB)) * ((Ea)/((1.+ma/mA)*Q+Ea)));
        double coscm=cos(cmangle);
        double coslab=(gamma+coscm)/sqrt(1+gamma*gamma+2*gamma*coscm);
        return acos(coslab)/d2r;

}

double Analysis::GetLabAngle16Op(double cmangle) {
        cmangle*=d2r;
        double ma=1;
        double mA=16.;
        double mb=4.;
        double mB=13;
        double Q=-2.2156-3.09;
        double Ea=Ebeam;
        double gamma=sqrt(((ma*mb)/(mA*mB)) * ((Ea)/((1.+ma/mA)*Q+Ea)));
        double coscm=cos(cmangle);
        double coslab=(gamma+coscm)/sqrt(1+gamma*gamma+2*gamma*coscm);
        return acos(coslab)/d2r;

}

double Analysis::GetCMAngle12C(double labangle) {
	double bestang=0;
	double bestdiff=180;
	for(double i=0;i<180;i+=0.1) {
		double thisang=GetLabAngle12C(i);
		if(fabs(thisang-labangle)<bestdiff) {
			bestdiff=fabs(thisang-labangle);
			bestang=i;
		}
	}
        return bestang;

}

double Analysis::GetCMAngle16O(double labangle) {
	double bestang=0;
	double bestdiff=180;
	for(double i=0;i<180;i+=0.1) {
		double thisang=GetLabAngle16O(i);
		if(fabs(thisang-labangle)<bestdiff) {
			bestdiff=fabs(thisang-labangle);
			bestang=i;
		}
	}
        return bestang;
}

double Analysis::GetCMAngle16Op(double labangle) {
	double bestang=0;
	double bestdiff=180;
	for(double i=0;i<180;i+=0.1) {
		double thisang=GetLabAngle16Op(i);
		if(fabs(thisang-labangle)<bestdiff) {
			bestdiff=fabs(thisang-labangle);
			bestang=i;
		}
	}
        return bestang;
}

