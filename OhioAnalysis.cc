//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Ohio Analysis:												//////
////		Can be executed either through root -l .X OhioAnalysis.cc+ or through the compiled version	//////
////		which uses int main()										//////
////		Set readmc to true to run MC data								//////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "OhioAnalysis.h"
using namespace std;
const double d2r=atan(1.)/45.;
const int npoints_dedx=500;
bool do_once=true;
bool server_mode=true;
int hitmult;
int main(int argc, const char * argv[]) {
	bool ismc=false;
	bool istwopmode=false;
	int experiment=1;
	server_mode=false;
	Analysis *hoyle = new Analysis(ismc);//Define new class object
	hoyle->DefineHists();//Define histograms
	hoyle->ReadMap();//Read the pixel to ASAD/AGET/CHN map from file
	hoyle->SetResponse2();//Read response file
	hoyle->SetExperiment(experiment);
	cout<<BOLDMAGENTA<<"Analysing Ohio data"<<RESET<<endl;
	TString arg1(argv[1]);
	TString arg2(argv[2]);
	int rnumb=atoi(argv[3]);
	cout<<"Run number "<<rnumb<<endl;
	if(!ismc)arg2+="_out.root";
	if(ismc) arg2+="_mc.root";
	cout<<arg1<<"\t"<<arg2<<endl;
	hoyle->SetHistName(arg2);
	hoyle->SetRunNumber(rnumb);
	hoyle->ReadCalibrationFile("alpha_calibration_253.txt");
	hoyle->Execute(arg1);//Option to read a single file
	cout<<hoyle->treetrack<<endl;
	hoyle->filetrack->cd();
	hoyle->treetrack->Write();
	cout<<hoyle->filetrack<<endl;
	hoyle->filetrack->Close();
	arg2="analysis_"+arg2;
	hoyle->WriteHist(arg2);//write the histogram name (want hist_{0}.root or something)
	return 0;
}

void OhioAnalysis() {//ROOT mode
//	Draw starting logo
        ifstream infile;
        infile.open("banner.txt");
	while(infile.good()) {
		string line;
		infile>>line;
		cout<<BOLDYELLOW<<line.c_str()<<RESET<<"\n";
	}
//	TCanvas *c1 = new TCanvas("c1","Canvas",0,0,0,0);
//	Define Canvas (0 size) for the TSpectrum which is called later
	bool ismc=false;
	bool istwopmode=false;
	int experiment=1;//
	server_mode=false;
	Analysis *hoyle = new Analysis(ismc);//Define new class object
	if(!ismc)hoyle->SetHistName("histo_Ohio.root");
	if(ismc) {
		if(experiment==1)hoyle->SetHistName("histo_mc_n2.root");
		if(experiment==2)hoyle->SetHistName("histo_mc_16Ona.root");
		if(experiment==3)hoyle->SetHistName("histo_mc_n0.root");
	}
	hoyle->DefineHists();//Define histograms
	hoyle->ReadMap();//Read the pixel to ASAD/AGET/CHN map from file
	hoyle->SetResponse2();
	hoyle->ReadCalibrationFile("alpha_calibration_253.txt");
	hoyle->SetExperiment(experiment);
	hoyle->twopmode=istwopmode;
	cout<<BOLDMAGENTA<<"Analysing Ohio data"<<RESET<<endl;
//	hoyle->Execute("/hdfs/user/jackbishop/run_0061.dat.12-03-19_07h48m24s.root");//
//	hoyle->Execute("online.root");//
	if(ismc) {
		if(experiment==1)hoyle->Execute("/hdfs/user/jackbishop/MM_nn2.root");
		if(experiment==2)hoyle->Execute("/hdfs/user/jackbishop/MM_16Ona.root");
		if(experiment==3)hoyle->Execute("/hdfs/user/jackbishop/MM_carbonel.root");
		if(experiment==3)hoyle->Execute("/hdfs/user/jackbishop/MM_oxygenel.root");
	}
	else {
		int errorcatch=0;
//		hoyle->Execute("/hdfs/user/jackbishop/run0320a/rootfiles/run_0223.dat.20-08-20_11h41m33s.root");
//		hoyle->Execute("/hdfs/user/jackbishop/run0320a/rootfiles/run_0228.dat.20-08-20_23h16m11s.root");
//		hoyle->Execute("/data/grgroup/gr02/jackbishop/run_0057.dat.10-03-20_15h36m23s.root");
//		hoyle->Execute("/data/grgroup/gr02/jackbishop/external/data/run0320a/acquisition/run/run_0096.dat.12-03-20_19h19m08s.13.root");
		errorcatch = hoyle->ExecuteList("Ed6500keV.txt");//Read a file with the name of files we want to analyse
//		errorcatch = hoyle->ExecuteList("run77.txt");//Read a file with the name of files we want to analyse
//		errorcatch = hoyle->ExecuteList("run81hdfs.txt");//Read a file with the name of files we want to analyse
//		errorcatch = hoyle->ExecuteList("batch_8.txt");//Read a file with the name of files we want to analyse
//		hoyle->Execute("/hdfs/user/jackbishop/run0320a/rootfiles/run_0203.dat.18-08-20_02h31m39s.root");
//		errorcatch = hoyle->ExecuteList("Ed5360keV_1.txt");//Read a file with the name of files we want to analyse
//		errorcatch = hoyle->ExecuteList("localfiles2.txt");//Read a file with the name of files we want to analyse
//		errorcatch = hoyle->ExecuteList("Ed5600keV.txt");//Read a file with the name of files we want to analyse
//		hoyle->Execute("/data/grgroup/gr02/jackbishop/MFMHistServer/run_0201.dat.17-08-20_22h25m48s.root");//Read a file with the name of files we want to analyse
//		hoyle->Execute("/data/grgroup/gr02/jackbishop/external/data/run0320a/acquisition/run/run_0042.dat.08-03-20_11h42m38s.root ");//Alpha test
//		hoyle->Execute("/data/grgroup/gr02/jackbishop/external/data/run0320a/acquisition/run/run_0040.dat.08-03-20_08h49m17s.root");//Alpha test
//		hoyle->Execute("/data/grgroup/gr02/jackbishop/external/data/run0320a/acquisition/run/run_0049.dat.10-03-20_07h13m05s.root");//Alpha test
//		hoyle->Execute("/data/grgroup/gr02/jackbishop/external/data/run0320a/acquisition/run/run_0041.dat.08-03-20_08h55m13s.root");//Alpha test
//		hoyle->Execute("/data/grgroup/gr02/jackbishop/MFMHistServer/run_0121.dat.06-08-20_12h43m35s.root");//pulser
//		hoyle->Execute("/hdfs/user/jackbishop/run0320a/rootfiles/run_0077.dat.10-03-20_17h57m37s.root");
//		errorcatch = hoyle->ExecuteList("alphacalib.txt");//Read a file with the name of files we want to analyse
//		errorcatch = hoyle->ExecuteList("runlist_good.txt");//Read a file with the name of files we want to analyse
		if(errorcatch==-1) {//If the list isn't found
			cout<<BOLDRED<<"An error occured! Check the run files"<<RESET<<endl;
			return;//Execute List failed
		}
	}
	cout<<GREEN<<"List completed\n"<<RESET;
//	TString outputdefaultname;
//	hoyle->Fits();
	hoyle->WriteHist(hoyle->GetHistName());//Write histograms to ROOT file
	cout<<BOLDGREEN<<"FINISHED\n"<<RESET;
	new TBrowser("historesponse.root");//Launch TBrowser when done
	return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis::Execute(TString fname) {
//	for(int i=0;i<180;i++) {
//		double CMtheta=i;
//		double labtheta=GetLabAngle12C(CMtheta);
//		double getback=GetCMAngle12C(labtheta);
//		cout<<CMtheta<<"\t"<<labtheta<<"\t"<<getback<<endl;
//	}

	if(Ebeam==0) Ebeam=9.5;
	cout<<"Ebeam = "<<Ebeam<<endl;
	DrawTheory();
	E_scaling = 0.03*(0.043/2400.);
	mc_thresh=0.;
	if(fexperiment==0) {
		cout<<BOLDRED<<"Don't forget to set the experiment to 1 for beta-decay or 2 for n,n'"<<RESET<<endl;
		return;
	}
	fexperiment=1;
	fulleventreadout=true;
	if(readmc || !twopmode) timetemp=512;
	timetemp=512;
	if(readmc) cout<<BOLDCYAN<<"MONTE CARLO MODE: TIMEBUCKETS = "<<RESET<<timetemp<<endl;
	debug=false;
	cout<<BOLDCYAN<<"Reading: "<<fname<<RESET<<endl;
//	Define parameters for reading of tree
	Int_t mmMul,evt,events,mmD2PTime;
	double answer,mean;
	TString ftrackname="tracks_"+GetHistName();
	filetrack = new TFile(ftrackname,"RECREATE");
	treetrack = new TTree("TATO","TATO track");
	treetrack->Branch("pointE",&pointE);
	treetrack->Branch("pointx",&pointx);
	treetrack->Branch("pointy",&pointy);
	treetrack->Branch("pointz",&pointz);
	treetrack->Branch("vertex",p_vertex,"vertex[3]/D");
	treetrack->Branch("end1",p_end1,"end1[3]/D");
	treetrack->Branch("end2",p_end2,"end2[3]/D");
	treetrack->Branch("chisqrd",&p_chisqrd);
	treetrack->Branch("En",&Ebeam);
	TFile *f;//Define file
	f=TFile::Open(fname);//Open ROOT file
	if(!f) {
		cout<<BOLDRED<<"COULD NOT OPEN FILE: "<<fname<<". Skipping..."<<RESET<<endl;
		return;
	}
        TFile *fcut;
        fcut = TFile::Open("maincut.root");
        gcut_1 = (TCutG*)fcut->Get("CUTG");

	const int timebuckets=timetemp;
	
/////////////// Read Tree/////////////////////////////////////////////////////////////////////////////////////////////////////////
	TTree *tree = new TTree;
	f->GetObject("TEvent",tree);//gets the tree from the ROOT file
//	tree->Scan("mmEventIdx:mmAsad:mmAget:mmChan:mmDecayNo:mmFrameNo");//Can print if more info is needed about the tree structure
	TBranch *bmult = 0,*bwavx=0,*bwavy=0,*bPixX=0,*bPixY=0,*bchan=0,*basad=0,*baget=0,*bcobo=0,*bdecay=0,*bD2P=0,*bhit=0,*bFrame=0,*bEvent=0;//define branches (not really needed)
//	Set branch addresses
//	tree->Print();
	tree->SetBranchAddress("mmMul",&mmMul,&bmult);
	tree->SetBranchAddress("mmHit",&mmHit,&bhit);
	tree->SetBranchAddress("mmCobo",&mmCobo,&bcobo);
	tree->SetBranchAddress("mmAsad",&mmAsad,&basad);
	tree->SetBranchAddress("mmAget",&mmAget,&baget);
	tree->SetBranchAddress("mmAsad",&mmAsad,&basad);
	tree->SetBranchAddress("mmChan",&mmChan,&bchan);
	tree->SetBranchAddress("mmDecayNo",&mmDecayNo,&bdecay);
	tree->SetBranchAddress("mmFrameNo",&mmFrameNo,&bFrame);
	tree->SetBranchAddress("mmEventIdx",&mmEventIdx,&bEvent);
	tree->SetBranchAddress("mmD2PTime",&mmD2PTime,&bD2P);
	tree->SetBranchAddress("mmWaveformY",&mmWaveformY);
	tree->SetBranchAddress("mmSCA1",&mmSCA1);
	tree->SetBranchAddress("mmSCA2",&mmSCA2);
	tree->SetBranchAddress("mmSCA3",&mmSCA3);
	tree->SetBranchAddress("mmSCA4",&mmSCA4);
	tree->SetBranchAddress("mmSCA5",&mmSCA5);
	tree->SetBranchAddress("mmSCA6",&mmSCA6);
	cout<<"Branches set"<<endl;
	events=tree->GetEntries();//load tree/
	cout<<"We have "<<BLUE<<events<<RESET<< " events\n";
	cout<<"Drawing "<<draw_events<<" events"<<endl;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	doublecount=0;
////////CORRECT MC MAPPING ERRORS FROM TEXATRESPONSE////////////////
	if(readmc && do_once) {
		for(int i=0;i<4;i++) {//ASAD
			for(int j=0;j<4;j++) {//AGET
				for(int k=0;k<64;k++) {//CHN
					if(map_x[i][j][k]>64 && map_x[i][j][k]<69 && ((map_y[i][j][k]>=48 && map_y[i][j][k]<128) || (map_y[i][j][k]>=16 && map_y[i][j][k]<32))) map_x[i][j][k]=133-map_x[i][j][k];
				}
			}
		}
		do_once=false;
	}
/////////////////////////////////////////
//	double responseshift[timebuckets];
//	int offset=90;
//	int source_offset=100;
	source_offset=100;
	if(readmc) source_offset=0;
        bool fpnadjust=false;
	if(fexperiment==1) fpnadjust=true;//
	if(fexperiment==2) fpnadjust=true;//
	if(fexperiment==3) fpnadjust=true;//

	bool takeall=true;
	if(takeall) cout<<BOLDRED<<"NOTE: Taking all events (no cuts applied)"<<RESET<<endl;

	int event_selection=0;
	if(verbose)cout<<"RESPONSE FUNCTION"<<endl;
	for(int t=0;t<timebuckets;t++) {
		responseshift[t]=response[(offset2+t+512)%timebuckets];
//		if(verbose)cout<<t<<"\t"<<response[t]<<endl;
	}
	outputdefaultname=GetHistName();
	starting_event=660;
////	events=101;
	const int limevents=500;
	if(events>limevents && false) {
		events=limevents;
		cout<<"Limiting events to "<<limevents<<endl;
	}
	rnd3 = new TRandom3();
	good=false;
	for(int i=starting_event;i<events;i++) {
//		if(i%10!=0) continue;
                pointx.clear();
                pointy.clear();
                pointz.clear();
                pointE.clear();

//		cout<<outputdefaultname<<endl;
//		cout<<"SCALER: "<<mmSCA1<<"\t"<<mmSCA2<<"\t"<<mmSCA3<<"\t"<<mmSCA4<<"\t"<<mmSCA5<<"\t"<<mmSCA6<<endl;
		if(debug) cout<<"Event: "<<i<<endl;
		if(i%50==0 && i>0) WriteHist(outputdefaultname);
		evtno=i;
		event_selection=0;
		totalEnergy=0;
		largest_time=512;
		hevents->Fill(event_selection);
		hevents->GetXaxis()->SetBinLabel(event_selection+1,"ALL");
		event_selection++;
		highy=0;
		lowy=128;
		if(good) {
			good_evt++;//last event was good so move on
			cout<<BOLDGREEN<<"Event "<<i-1<<" was good"<<RESET<<endl;
		}
		DefineEventHistograms();
		good=false;
		fulleventreadout=false;
		if(debug)cout<<"____________________________NEW EVENT__________________________________________________"<<endl;
		tree->GetEntry(i);//Get the entries from the tree
		htime->Fill(mmD2PTime);
		htime2->Fill(mmD2PTime);
		hitmult=mmHit;
		mmHit=mmMul;
		if(mmHit>1000) fulleventreadout=true;
		if(verbose)cout<<"Mult = "<<mmHit<<endl;
/*
		if(mmHit<=1088) {
			good=false;
			continue;
		}
*/
		if(fulleventreadout && i%100==0 && server_mode) cout<<"Event: "<<i<<" out of "<<events<<" in full-event mode readout\n";
		if(fulleventreadout && i%10==0 && !server_mode) cout<<"Event: "<<i<<" out of "<<events<<" in full-event mode readout\n";
		else{
			if(i%250==0)cout<<"Event: "<<i<<" out of "<<events<<"\n";
		}
//		if(mmD2PTime<200e3 && !readmc && twopmode) continue;//cut on d2p time - if time between triggers is too small then reject (probably noise and can be correct for after) 100us is actual value - cut above 200us
                hevents->Fill(event_selection);
                hevents->GetXaxis()->SetBinLabel(event_selection+1,"Passed D2P cut");
                event_selection++;
		if(debug)cout<<"Done"<<endl;
		if(good_evt==draw_events&& good_evt!=events-1) cout<<BLUE<<"Drawing no more good events ("<<(draw_events)<<")\n"<<RESET;
		double padenergy,peakEnergy=0,peaking,peaktime=0,sumE;
		if(debug) cout<<"Point Echo\n";
		const int mmhits=2400;//maximum MM hits to consider (arbitrarly high choice)
		double mmenergy[mmhits][2],mmtime[mmhits][2],mmenergy_2[mmhits][2];
		if(debug)cout<<"Big arrays"<<endl;
//		const int striparsize=140;//make this array large enough to fit the width of the pseudopixels (140)
//		double mmstripsE[striparsize][2][2],mmchainsE[striparsize][2],mmstripsT[striparsize][2][2],mmchainsT[striparsize][2],mmstripsE_2[striparsize][2][2],mmchainsE_2[striparsize][2];
		for(int x=0;x<striparsize;x++) {
			for(int y=0;y<2;y++) {//zero arrays
				for(int k=0;k<2;k++) {
					mmchainsE[x][y]=0.;
					mmstripsE[x][k][y]=0.;
					mmchainsT[x][y]=0.;
					mmstripsT[x][k][y]=0.;
					mmstripsE_2[x][k][y]=0.;
					mmchainsE_2[x][y]=0.;
//					cout<<x<<"\t"<<y<<"\t"<<k<<"\t"<<striparsize<<endl;
				}
			}
		}
		if(debug)cout<<"Point Foxtrot\n";
		if(hcurrenttrack3D!=0) delete hcurrenttrack3D;//If hcurrenttrack3D is still hanging around, get rid of it as we want to make a new one
		if(hcurrenttrack3D_2!=0) delete hcurrenttrack3D_2;//If hcurrenttrack3D_2 is still hanging around, get rid of it as we want to make a new one
		if(hcurrenttrack3D_3!=0) delete hcurrenttrack3D_3;//If hcurrenttrack3D_3 is still hanging around, get rid of it as we want to make a new one
		if(hcurdedx!=0) delete hcurdedx;
		if(hangle!=0) delete hangle;
		if(debug)cout<<"point Foxtrot 2\n";
		if(hbeam_c!=0) delete hbeam_c;
		if(hflat_all!=0) delete hflat_all;
		if(hHoughXYb!=0) delete hHoughXYb;
		if(hHoughXZb!=0) delete hHoughXZb;
		if(hHoughYZb!=0) delete hHoughYZb;
		if(hcurrentYZ!=0) delete hcurrentYZ;
		hcurrenttrack3D = new TH3F("hcurrenttrack3D","Current track",140,0,140,128,0,128,timebuckets,0,timebuckets);//3D track
		hcurrenttrack3D_2 = new TH3F("hcurrenttrack3D_2","Current track - peak energies",140,0,140,128,0,128,timebuckets,0,timebuckets);//3D track
		hcurrenttrack3D_3 = new TH3F("hcurrenttrack3D_3","Current track long arm",140,0,140,128,0,128,timebuckets,0,timebuckets);//3D track
		hflat_all = new TH2F("hflat_all","Flattened distance",100,-100,100,100,-100,100);//flattened track
		hangle = new TH1F("hangle","Dot product for Ys",45,0,180);
		hcurdedx = new TH1F("hcurdedx","Current dedx",npoints_dedx,0.,250.);
		hbeam_c = new TH2F("hbeam_c","BEAM",140,0,140,128,0,128);

                hHoughXYb = new TH2F("hHoughXYb","Beam Hough XY",1800,0,180,600,-50,50);
                hHoughXZb = new TH2F("hHoughXZb","Beam Hough XZ",1800,0,180,600,-50,50);
                hHoughYZb = new TH2F("hHoughYZb","Beam Hough YZ",1800,0,180,600,-50,50);

		hcurrentYZ = new TH2F("hcurrentYZ","",128,0,128,timebuckets,0,timebuckets);
		double Ebeam=0;
		double firstpadE[6]={0.};
		double secondpadE[6]={0.};
		int mult_decay=0;
///////////////////////////////////////////////////////////////////////////////////////////////////////////
		int firstE=0,secondE=0,second_E=0;
		bool mixedevent=false;
		int largesty=0,largesty_b=0;
//////////////////////////////////////////////////////////////////////////////////////////////////////////
////		Proper analysis starts here								//
//////////////////////////////////////////////////////////////////////////////////////////////////////////
////		Analyze the waveforms to get energy and time information				///
///////////////////////////////////////////////////////////////////////////////////////////////////////////
		totalEnergy_glob=0;
		if(verbose || debug)cout<<"Event "<<i<<" or "<<mmEventIdx<<"  has "<<mmHit<<" "<<mmMul<<" mmHits"<<endl;
		bool escape=false;
		bool reject_escape=false;
		int frontsimult=0,backsimult=0,frontsichan=0,backsichan=0;
		double frontsiE=0.,backsiE=0.;
	        double fpnsum[512][2]={0.},fpnaverage[2]={0.},dEleft=0.,dEright=0.;
	        for(int j=0;j<mmHit;j++) {//find out what those pesky fpns are up to
	                if(mmCobo[j]==1 && (mmChan[j]==11 || mmChan[j]==22 ||mmChan[j]==45 ||mmChan[j]==56)) {
	                        for(int t=0;t<512;t++) {
	                                fpnsum[t][mmAget[j]]+=mmWaveformY[j][t];
	                                fpnaverage[mmAget[j]]+=mmWaveformY[j][t];
	                        }
			}
	        }
		step_count=0;
		double largest_width=0.;
		if(good_evt<draw_events) hpulse[good_evt]->SetTitle(Form("Mult: %d",mmHit));
		if(verbose || debug) cout<<"Looping over "<<mmHit<<" waveforms"<<endl;
		for(int j=0;j<mmHit;j++) {//loop over waveforms
				waveformno=j;
				if(mmChan[j]==11 || mmChan[j]==22 || mmChan[j]==45 || mmChan[j]==56) continue;
	                        if(mmCobo[j]==1)for(int t=0;t<512;t++) mmWaveformY[j][t]+=-fpnsum[t][mmAget[j]]*0.25+((fpnaverage[mmAget[j]]*0.25)/512.);//subtract FPN baseline

	                        if(fpnadjust) {//adjust channels to account for FPN
        	                        double chan,dchan=0;
                	                chan=mmChan[j];
                        	        if(chan<11) {
                               	        	dchan = chan;
                               		}else if(chan>11 && chan<22) {
                                	        dchan = chan - 1;
                                	}else if(chan>22 && chan<45) {
                                	        dchan = chan - 2;
                                	}else if(chan>45 && chan<56) {
                                	        dchan = chan - 3;
                                	}else if(chan>56) {
                                	        dchan = chan - 4;
                                	}
                                	mmChan[j]=dchan;//new channel
                        	}
				mmDecayNo[j]=1;
//				hwaveform[j]->Reset();
				if(mmHit==1098)AnalyzeSi();
				double maxE=RemoveBG();
				for(int t=0;t<timebuckets;t++) {
					if(good_evt<draw_events && maxE>50) hpulse[good_evt]->Fill(t,mmWaveformY[j][t]);
				}
                                double maxT=0;
	                        if(mmCobo[j]==0 && (map_y[mmAsad[j]][mmAget[j]][mmChan[j]]==-1 || map_x[mmAsad[j]][mmAget[j]][mmChan[j]]==-1) && good_evt<draw_events) {
        	                        int side = (mmAsad[j]==3) ? 0 : 1;//if ASAD==3 then left side of beam
                	                int lowx,highx;
                        	        lowx=(side==1)?0:68;
                                	highx=(side==1)?64:128;
					if(side==0) dEleft+=maxE;
					if(side==1) dEright+=maxE;
                                	for(int xpos=lowx;xpos<highx;xpos++) {
                                	        hmaxstrips[good_evt]->Fill(xpos,map_y[mmAsad[j]][mmAget[j]][mmChan[j]],maxE);
                                	}
                                	for(int ypos=0;ypos<140;ypos++) {
                                	        hmaxchains[good_evt]->Fill(map_x[mmAsad[j]][mmAget[j]][mmChan[j]],ypos,maxE);
                                	}
                        	}

	                        if(mmCobo[j]==0 && (map_y[mmAsad[j]][mmAget[j]][mmChan[j]]!=-1 && map_x[mmAsad[j]][mmAget[j]][mmChan[j]]!=-1) && good_evt<draw_events) {
                                        hmaxchains[good_evt]->Fill(map_x[mmAsad[j]][mmAget[j]][mmChan[j]],map_y[mmAsad[j]][mmAget[j]][mmChan[j]],maxE);
        	                }
//				cout<<j<<"\t"<<maxE<<endl;
				if(maxE<150) {
					delete hwave;
					delete hwave2;
					continue;
				}
				if(debug) cout<<"E: "<<maxE<<endl;
				totalEnergy+=integratewf;
//				cout<<"MAXE = "<<maxE<<endl;
				double tempmax=0;
				for(int t=0;t<timebuckets;t++) {
					if(good_evt<draw_events && maxE>50) hpulse_b[good_evt]->Fill(t,mmWaveformY[j][t]);
					if(mmWaveformY[j][t]>tempmax) {
						tempmax=mmWaveformY[j][t];
						maxT=t;
					}
				}

				//plot silicon energies
//				if(mmCobo[j]==1 && maxE>200 && ((mmChan[j]>=32 && mmChan[j]<49 && mmAget[j]==0) || (mmChan[j]>=32 && mmChan[j]<36 && mmAget[j]==1))) {
				if(mmCobo[j]==1 && maxE>200 && ((mmChan[j]>=0 && mmChan[j]<68 && mmAget[j]==0) || (mmChan[j]>=0 && mmChan[j]<68 && mmAget[j]==1))) {
					int silicon_channel=0;
					silicon_channel=mmChan[j]-32+mmAget[j]*16;//0->15 are fronts then 16,17,18,19 are backs
//					silicon_channel=mmChan[j]+mmAget[j]*64;//0->15 are fronts then 16,17,18,19 are backs
//					cout<<"SILICON HIT "<<silicon_channel<<endl;
					hSilicon[silicon_channel]->Fill(maxE);
					if(mmAget[j]==0) {
						frontsimult++;
						frontsiE=maxE;
						frontsichan=silicon_channel;
//						cout<<"Front si "<<mmChan[j]<<endl;
					}
					if(mmAget[j]==1) {
						backsimult++;
						backsiE=maxE;
						backsichan=silicon_channel;
//						cout<<"Back si "<<mmChan[j]<<endl;
					}
				}
				peakEnergy+=maxE;//add total energy
				hallhits->Fill(map_x[mmAsad[j]][mmAget[j]][mmChan[j]],map_y[mmAsad[j]][mmAget[j]][mmChan[j]]);
				if(mmCobo[j]==1) {
					delete hwave;
					delete hwave2;
					continue;//Silicons are done here
				}
				if(mmCobo[j]==0 && maxE>150) mult_decay++;
				if(reject_escape && ((map_y[mmAsad[j]][mmAget[j]][mmChan[j]]==0 || map_y[mmAsad[j]][mmAget[j]][mmChan[j]]==127) || map_y[mmAsad[j]][mmAget[j]][mmChan[j]]==1 || map_y[mmAsad[j]][mmAget[j]][mmChan[j]]==126)) {
//					good=false;
					if(verbose)cout<<"Escaped"<<endl;
					escape=true;
				}
				if(map_x[mmAsad[j]][mmAget[j]][mmChan[j]]!=-1 && map_y[mmAsad[j]][mmAget[j]][mmChan[j]]!=-1) {//fill central region
					centralhitsevt[map_x[mmAsad[j]][mmAget[j]][mmChan[j]]-64][map_y[mmAsad[j]][mmAget[j]][mmChan[j]]]=maxE;
					centralhitsTevt[map_x[mmAsad[j]][mmAget[j]][mmChan[j]]-64][map_y[mmAsad[j]][mmAget[j]][mmChan[j]]]=maxT;
				}
				if(good_evt<draw_events){//draw raw hits (no cleaning)
					hHits[good_evt]->Fill(map_x[mmAsad[j]][mmAget[j]][mmChan[j]],map_y[mmAsad[j]][mmAget[j]][mmChan[j]],maxE);
					hHits_side[good_evt]->Fill(map_y[mmAsad[j]][mmAget[j]][mmChan[j]],maxT,maxE);
				}
				if(map_x[mmAsad[j]][mmAget[j]][mmChan[j]]==-1) {//strips
					if(mmAsad[j]==2)for(int l=0;l<63;l++) {
						hallhits->Fill(l,map_y[mmAsad[j]][mmAget[j]][mmChan[j]]);
						if(good_evt<draw_events) hHits[good_evt]->Fill(l,map_y[mmAsad[j]][mmAget[j]][mmChan[j]],maxE);
					}
					if(mmAsad[j]==3)for(int l=70;l<140;l++) {
						hallhits->Fill(l,map_y[mmAsad[j]][mmAget[j]][mmChan[j]]);
						if(good_evt<draw_events) hHits[good_evt]->Fill(l,map_y[mmAsad[j]][mmAget[j]][mmChan[j]],maxE);
					}
				}
				if(map_y[mmAsad[j]][mmAget[j]][mmChan[j]]==-1) {//chains
					for(int l=0;l<128;l++) {
						hallhits->Fill(map_x[mmAsad[j]][mmAget[j]][mmChan[j]],l);
						if(good_evt<draw_events) hHits[good_evt]->Fill(map_x[mmAsad[j]][mmAget[j]][mmChan[j]],l,maxE);
					}
				}
				peaking=0;
				sumE=0;
				if(debug)cout<<"Point Golf\n";
				bool simple=false;//change to true to use peak finding method
				for(int m=0;m<10;m++) xgoodpeaks[m]=0.;
				good_double=0;
                                if(simple || saturated) {
					good_double=1;
                                        double maxT2=0.,maxE2=0.;
                                        for(int t=0;t<timebuckets-100;t++) {
                                                if(hwave->GetBinContent(t)>maxE2) {
                                                        maxE2=hwave->GetBinContent(t);
                                                        maxT2=t-source_offset;
                                               }
                                        }
                                        xgoodpeaks[0]=maxT2;
					if(debug) cout<<"Peaking at t= "<<maxT2<<endl;
//					if(saturated) cout<<"Saturated, peaking at t= "<<maxT<<" with "<<maxE<<endl;
					if(saturated) {
						FitSaturated(maxE,maxT,maxT2,hwave2);
					}
                                }
				else{
					if(debug)cout<<"CONV"<<endl;
					Deconvolution();
				}
				if(debug)cout<<"Done deconv"<<endl;
//				Fit to see if multiple tracks here from width of deconvoluted spectrum
				for(int l=0;l<good_double;l++) {
					mygaussian->SetParameter(0,hwave2->GetBinContent(xgoodpeaks[l]+1));
					mygaussian->SetParameter(1,xgoodpeaks[l]);
					mygaussian->SetParameter(2,10.);
					mygaussian->SetParameter(3,0.);
					mygaussian->SetParameter(4,0.);

					int statusfit = hwave2->Fit(mygaussian,"Q");
					if(fabs(mygaussian->GetParameter(1)-xgoodpeaks[l])<10) {
						xgoodpeaks[l]=mygaussian->GetParameter(1);
						double widthgauss=mygaussian->GetParameter(2);
						if(widthgauss>largest_width) largest_width=widthgauss;
					}
				}
//
				hpeaks->Fill(good_double);
				if(debug)cout<<"sorted"<<endl;
				if(debug)cout<<"energy,time for "<<good_double<<" peaks"<<endl;
				double sep=50;
				if(good_double==2)sep=xgoodpeaks[1]-xgoodpeaks[0];
				if(good_double<=1) {//If zero or one hit
						double maxE2=maxE;
						double maxT2=maxT;
						if(good_double==1) {
							maxE2=mmWaveformY[j][int(xgoodpeaks[0])];
							maxT2=xgoodpeaks[0];
						}
						FitSaturated(maxE,maxT,maxT2,hwave2);
						mmenergy_2[j][0]=maxE;
						mmenergy[j][0]=maxE;
						mmtime[j][0]=maxT;
						good_double=1;//we now have one good peak
						xgoodpeaks[0]=maxT;
					}

				for(int l=0;l<min(2,good_double);l++) {//allocate/draw multiple peaks for each channel
					if(simple)mmenergy[j][l]=hwave->GetBinContent(xgoodpeaks[l]+1);
//					cout<<mmWaveformY[j][int(xgoodpeaks[l]-source_offset)]<<endl;
//					mmenergy_2[j][l]=mmWaveformY[j][int(xgoodpeaks[l]-source_offset)];//take the peak value from the raw waveform
					mmenergy_2[j][l]=mmWaveformY[j][int(xgoodpeaks[l])];//take the peak value from the raw waveform
					mmenergy[j][l]=mmenergy_2[j][l];
					if(mmenergy[j][l]<70) continue;
					mmtime[j][l]=int(xgoodpeaks[l]);
//					cout<<"Time = "<<mmtime[j][l]<<"\t"<<mmenergy[j][l]<<endl;
					if(debug) cout<<"Normal energy: "<<mmenergy[j][l]<<" for "<<xgoodpeaks[l]<<endl;
					if(debug) cout<<l<<"\tof\t"<<good_double<<"\t"<<xgoodpeaks[l]<<endl;
					if(good_double==1) {
						double maxE2=mmenergy_2[j][0];
						double maxT2=mmtime[j][0];
						FitSaturated(maxE,maxT,maxT2,hwave2);
						mmenergy_2[j][l]=maxE;
						mmenergy[j][l]=maxE;
						mmtime[j][l]=maxT;
					}
					if(good_evt<draw_events) {
//						cout<<"Waveform hit"<<good_evt<<"\t"<<draw_events<<"\t"<<xgoodpeaks[l]-source_offset<<"\t"<<mmWaveformY[j][int(xgoodpeaks[l])]<<"\t"<<j<<endl;
						hwaveform_hit[good_evt]->Fill(j,mmtime[j][l],mmenergy_2[j][l]);
					}
					if(map_y[mmAsad[j]][mmAget[j]][mmChan[j]]>highy && mmenergy[j][l]>100.) highy=map_y[mmAsad[j]][mmAget[j]][mmChan[j]];
					if(map_y[mmAsad[j]][mmAget[j]][mmChan[j]]<lowy && map_y[mmAsad[j]][mmAget[j]][mmChan[j]]!=-1 && mmenergy[j][l]>100.) lowy=map_y[mmAsad[j]][mmAget[j]][mmChan[j]];
					if((mmtime[j][l]<=100 || (map_x[mmAsad[j]][mmAget[j]][mmChan[j]]==0) || (map_x[mmAsad[j]][mmAget[j]][mmChan[j]]==139) ||
						 (map_y[mmAsad[j]][mmAget[j]][mmChan[j]]==0) || (map_y[mmAsad[j]][mmAget[j]][mmChan[j]]==127)) && !takeall) good=false;//Event is escaping so we must discard it
//					if(mmtime[j][l]>300)cout<<saturated<<"\t"<<l<<"\t"<<good_double<<"\t"<<xgoodpeaks[l]<<"\t"<<mmenergy[j][l]<<"\t"<<mmtime[j][l]<<endl;
					if(good_evt<draw_events && map_x[mmAsad[j]][mmAget[j]][mmChan[j]]!=-1 && map_y[mmAsad[j]][mmAget[j]][mmChan[j]]!=-1) {//central pads
//						htrack3D[good_evt]->Fill(map_x[mmAsad[j]][mmAget[j]][mmChan[j]],map_y[mmAsad[j]][mmAget[j]][mmChan[j]],mmtime[j][l],mmenergy[j][l]);
						htrackxy[good_evt]->Fill(map_x[mmAsad[j]][mmAget[j]][mmChan[j]],map_y[mmAsad[j]][mmAget[j]][mmChan[j]],mmenergy[j][l]);
						hsuper[good_evt]->Fill(map_x[mmAsad[j]][mmAget[j]][mmChan[j]],map_y[mmAsad[j]][mmAget[j]][mmChan[j]]);
						htrackyz[good_evt]->Fill(map_y[mmAsad[j]][mmAget[j]][mmChan[j]],mmtime[j][l],mmenergy[j][l]);
						htrackxz[good_evt]->Fill(map_x[mmAsad[j]][mmAget[j]][mmChan[j]],mmtime[j][l],mmenergy[j][l]);
					}
					if(map_x[mmAsad[j]][mmAget[j]][mmChan[j]]!=-1 && map_y[mmAsad[j]][mmAget[j]][mmChan[j]]!=-1) {//central pads
						hcurrenttrack3D->Fill(map_x[mmAsad[j]][mmAget[j]][mmChan[j]],map_y[mmAsad[j]][mmAget[j]][mmChan[j]],mmtime[j][l],mmenergy[j][l]);//plot central pads as 3D track
						hcurrenttrack3D_2->Fill(map_x[mmAsad[j]][mmAget[j]][mmChan[j]],map_y[mmAsad[j]][mmAget[j]][mmChan[j]],mmtime[j][l],mmenergy_2[j][l]);//plot central pads as 3D track
						hcurrentYZ->Fill(map_y[mmAsad[j]][mmAget[j]][mmChan[j]],mmtime[j][l],mmenergy_2[j][l]);
						if(good_evt<draw_events) htrackxysmooth[good_evt]->Fill(map_x[mmAsad[j]][mmAget[j]][mmChan[j]],map_y[mmAsad[j]][mmAget[j]][mmChan[j]],mmenergy_2[j][l]);//weight by energy
					}
					if(map_x[mmAsad[j]][mmAget[j]][mmChan[j]]==-1) {//chains
						int side = (mmAsad[j]==3) ? 1 : 0;//if ASAD==3 then left side of beam
						if(good_evt<draw_events) htrackxz[good_evt]->Fill(map_x[mmAsad[j]][mmAget[j]][mmChan[j]],mmtime[j][l],mmenergy[j][l]);
						mmstripsE[map_y[mmAsad[j]][mmAget[j]][mmChan[j]]][side][l]=mmenergy[j][l];
						mmstripsE_2[map_y[mmAsad[j]][mmAget[j]][mmChan[j]]][side][l]=mmenergy[j][l];
						mmstripsT[map_y[mmAsad[j]][mmAget[j]][mmChan[j]]][side][l]=xgoodpeaks[l];
//						cout<<"CHAINS: "<<i<<"\t"<<j<<"\t"<<mmenergy[j][l]<<"\t"<<maxT<<endl;

					}
					if(map_y[mmAsad[j]][mmAget[j]][mmChan[j]]==-1) {//strips
						if(good_evt<draw_events)htrackyz[good_evt]->Fill(map_y[mmAsad[j]][mmAget[j]][mmChan[j]],mmtime[j][l],mmenergy[j][l]);
						mmchainsE[map_x[mmAsad[j]][mmAget[j]][mmChan[j]]][l]=mmenergy[j][l];
						mmchainsE_2[map_x[mmAsad[j]][mmAget[j]][mmChan[j]]][l]=mmenergy[j][l];
						mmchainsT[map_x[mmAsad[j]][mmAget[j]][mmChan[j]]][l]=xgoodpeaks[l];
//						cout<<"STRIPS: "<<i<<"\t"<<j<<"\t"<<mmenergy[j][l]<<"\t"<<maxT<<endl;

					}
					peaking+=hwave->GetBinContent(xgoodpeaks[l]);
				}
				if(debug) cout<<"beep boop"<<endl;
				delete hwave;//done with this so delete it
				delete hwave2;
//				if(step_count>500) {//This is a jumping waveform event - not real
//					good=false;
//					if(verbose) cout<<"Jump waveform "<<i<<endl;
//					j=mmHit;
//					continue;
//				}

		}//finish loop over all hits
		if(step_count>500) {//This is a jumping waveform event - not real
			if(verbose)cout<<"Jumping waveform "<<i<<" mult = "<<step_count<<endl;
			good=false;
			continue;
		}
//		cout<<"MULTS "<<frontsimult<<"\t"<<backsimult<<endl;
		if(debug) cout<<"Done loop"<<endl;
                hevents->Fill(event_selection);
                hevents->GetXaxis()->SetBinLabel(event_selection+1,"Jump check");
                event_selection++;

		hwidth->Fill(largest_width);

                h1->Fill(totalEnergy);//integrated energy
                h2->Fill(peakEnergy);//energy in keV
		hE1E2->Fill(totalEnergy,peakEnergy);
                hevents->Fill(event_selection);
                hevents->GetXaxis()->SetBinLabel(event_selection+1,"Good waveform - no jump");
                event_selection++;
		if(frontsimult==1 && backsimult==1) {
			hSilicong[frontsichan]->Fill(frontsiE);
			hSilicong[backsichan]->Fill(backsiE);
//			cout<<"Coincidence "<<frontsichan<<"\t"<<backsichan<<endl;
//			if(frontsiE>5 && backsiE>5) good=true;
			double bigger=0;
			if(dEleft>dEright) {
				bigger=dEleft;
			}
			else {
				bigger=dEright;
			}
			hdEE->Fill(bigger,frontsiE);
			if(mult_decay<4) {
				hSilicongg[frontsichan]->Fill(frontsiE*gain[frontsichan]+offset[frontsichan]);
				hSilicongg[backsichan]->Fill(backsiE*gain[backsichan]+offset[backsichan]);
//				cout<<"No MM"<<endl;
			}
		}
//		cout<<i<<"\t"<<mult_decay<<" decays"<<endl;

		if(hcurrenttrack3D->GetEntries()<5) {//too few
			good=false;
			if(verbose)cout<<"MM mult<5 for event "<<i<<"("<<hcurrenttrack3D->GetEntries()<<")"<<endl;
			continue;
		}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///		Use the energy and times to build strip/chain information				///
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
		fexperiment=2;//default to n,n0
//		cout<<evtno<<"\t"<<peakEnergy<<"\t"<<mult_decay<<endl;
		if(mult_decay>5) {
			FitLine();
		}
		if(pressure==100 && mult_decay<15) {
			good=false;
//			good=true;
//			Ransac2();
			continue;
		}
		good=true;
		if(mult_decay<30 && pressure==50) {//only select multitrack events
			good=false;
			fexperiment=3;//n,n0
			continue;
		}
		else {
			good=true;
//			fexperiment=2;//n,a or Hoyle
			fexperiment=1;//Hoyle mode
			h3->Fill(totalEnergy);
//			good=true;
		}
//		if(pointE.size()<30) {
//			good=false;
//			continue;
//		}
		cout<<"Event: "<<evtno<<" MULT: "<<mult_decay<<endl;

		if(debug) cout<<"Point here"<<endl;
		totalEnergy_glob=totalEnergy;
		hmult->Fill(mult_decay);
		spread=highy-lowy;
		if(spread<0) spread=0;
//		cout<<highy<<"\t"<<lowy<<endl;
		hmultspread->Fill(mult_decay,spread);
		hmultE->Fill(mult_decay,peakEnergy);
		ingate=gcut_1->IsInside(mult_decay,peakEnergy);
		if(verbose) cout<<mult_decay<<"\t"<<spread<<endl;
/*
		if(fexperiment==2) {
			if(mult_decay<80 && spread<35) {
				fexperiment=1;//Probably Hoyle
	//			Check if Hoyle
			}
			else {
				fexperiment=2;//default guess to (n,a)
			}
		}
*/
		if(mult_decay>200) {
			if(verbose) cout<<"Mult too high - "<<mult_decay<<endl;
			good=false;
			continue;
		}
//		cout<<peakEnergy<<endl;
//		fexperiment=3;//n,n0 default
		if(good_evt<draw_events) htrackxy[good_evt]->SetTitle(Form("Mult spread=%d %g",mult_decay,spread));
		if(good_evt<draw_events) htrackyz[good_evt]->SetTitle(Form("Mult E%d %g",mult_decay,peakEnergy));

		if(good_evt<draw_events)htrackxysmooth[good_evt]->SetTitle(Form("Max/min=%g %g",highy,lowy));//record the total energy on the 3D track histogram

		if(!good && !takeall && !readmc) continue;
		if(reject_escape && escape) {
			cout<<"Escaped "<<i<<endl;
			good=false;
			continue;
		}
                hevents->Fill(event_selection);
                hevents->GetXaxis()->SetBinLabel(event_selection+1,"HAS DECAY & DOESN'T ESCAPE");
                event_selection++;

		if(debug) cout<<"Clean strips/chains"<<endl;
		StripsChains();//Match strips and chains
                double dist_diff=0;
                int beamx=0,beamy=0,beamz=0;
		bool anytrack=false;
		bool goodarms=false;
		int missing_strip=0;
		if(debug) cout<<"Fill beam origin\t"<<beamx<<"\t"<<beamy<<endl;
                hevents->Fill(event_selection);
                hevents->GetXaxis()->SetBinLabel(event_selection+1,"GOOD BEAM");
                event_selection++;
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///		Clean up by using nearest method cut (<20) and plot the total energies
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////Nearest neighbour//////////////////////////////////////////////
//		NearestNeighbour();//Remove outliers
//		if(debug) cout<<"Outliers removed"<<endl;
		if(totalEnergy<0) cout<<"-ve energy:\t"<<i<<endl;
                hEx_evt->Fill(totalEnergy*E_scaling+7.272,100.*i/events);
                hEx1Ex2->Fill(totalEnergy*E_scaling+7.272,peakEnergy*(.043/2400.)+7.272);
                hEtotntot->Fill(totalEnergy,mmHit);
                hpeaktotal->Fill(peakEnergy,totalEnergy);
                hevents->Fill(event_selection);
                hevents->GetXaxis()->SetBinLabel(event_selection+1,"INSIDE ENERGY GATE");
                event_selection++;

                if(hcurrenttrack3D->GetEntries()==0 && !takeall) continue;
                hevents->Fill(event_selection);
                hevents->GetXaxis()->SetBinLabel(event_selection+1,"HAS GOOD 3D TRACK");
                event_selection++;

////////////////////////////////////////////////////////////
//		Find the center of the decay
////////////////////////////////////////////////////////////
		if(debug)cout<<"Find start/centre"<<endl;
		int maxvalue=0.,maxx=-1,maxy=-1,maxz=-1;
		GetEnd(dmaxx,dmaxy,dmaxz);//max(xyz) is center, 4. is diffuseness of Gaussian
		maxx=dmaxx;
		maxy=dmaxy;
		maxz=dmaxz;
		RemoveOutliers(maxx,maxy,maxz);
		if(maxy<5) {
			good=false;
			if(verbose) cout<<"Too close to start of MM "<<maxy<<" - jump mult = "<<step_count<<endl;
			continue;
		}
		Dist(maxx,maxy,maxz,centr_dist_x,centr_dist_y,centr_dist_z);
		Dist(maxx,maxy,maxz,centr_dist_x2,centr_dist_y2,centr_dist_z2);
//		cout<<"Start :"<<centr_dist_x<<"\t"<<centr_dist_y<<"\t"<<centr_dist_z<<endl;
		if(verbose) {
			cout<<"Classification: ";
			if(fexperiment==1) cout<<"(n,n2)"<<endl;
			if(fexperiment==2) cout<<"(n,a)"<<endl;
			if(fexperiment==3) cout<<"(n,n0)"<<endl;
		}

		if(fexperiment==2) {
			double highestE=0;
			for(int x=0;x<140;x++) {
				for(int y=0;y<128;y++) {
					for(int z=0;z<timetemp;z++) {
						if(hcurrenttrack3D_2->GetBinContent(x+1,y+1,z+1)>highestE) {
							highestE=hcurrenttrack3D_2->GetBinContent(x+1,y+1,z+1);
							maxx=x;
							maxy=y;
							maxz=z;
						}
					}
				}
			}
//			Find the interaction point, highest E dep
//			dmaxx=maxx;
//			dmaxy=maxy;
//			dmaxz=maxz;
		}

		if(fexperiment==2) {//16O(n,a)13C
			Ransac2();
//			Check if just a straight line first
			double chisqrdline=FitLine();
//			cout<<"chisqrd line= "<<chisqrdline<<endl;
			if(good_evt<draw_events)hdistxz1[good_evt]->SetTitle(Form("chisqrd %g",chisqrdline));
			if(chisqrdline<2) {//just a single track
				good=false;
				continue;
			}
			cout<<"event: "<<evtno<<endl;
			double neutronE12C=0.;
			double neutronE16O=0.;
			double chisqrd12C=1e3,chisqrd16O=1e3,chisqrd16Op=1e3;
			double chisqrd212C=1e3,chisqrd216O=1e3,chisqrd216Op=1e3;
			select_nalpha=1;
			chisqrd12C=Fitnalpha();
			chisqrd212C=chisqrd_nalpha;
			hchisqrdth12C->Fill(chisqrd_nalpha_theta);
			hchisqrdE12C->Fill(chisqrd_nalpha_dist);
			select_nalpha=2;
			chisqrd16O=Fitnalpha();
			chisqrd216O=chisqrd_nalpha;
			hchisqrdth16O->Fill(chisqrd_nalpha_theta);
			hchisqrdE16O->Fill(chisqrd_nalpha_dist);

			select_nalpha=3;
			chisqrd16Op=Fitnalpha();
			chisqrd216Op=chisqrd_nalpha;
			hchisqrdth16Op->Fill(chisqrd_nalpha_theta);
			hchisqrdE16Op->Fill(chisqrd_nalpha_dist);
//			if(pressure==100) chisqrd16Op=1e3;//Do not allow for 16O for 100 Torr measurements - cannot differentiate well
			double chisqrdcut=2.;

			hchisqrd->Fill(chisqrd12C,chisqrd16O);
			hselect2->Fill((chisqrd12C-chisqrd16O)/(chisqrd12C+chisqrd16O),min(chisqrd12C,chisqrd16O));
			if(good_evt<draw_events)hdistxz2a[good_evt]->SetTitle(Form("chisqrd (12C/16O/16O*) %g %g %g",chisqrd12C,chisqrd16O,chisqrd16Op));
//			chisqrd12C=dtheta[0];
//			chisqrd16O=dtheta[1];
//			cout<<fit_theta[0][0]<<"\t"<<fit_theta[0][1]<<"\t\t"<<fit_theta[1][0]<<"\t"<<fit_theta[1][1]<<endl;
			if(good_evt<draw_events)hdistxy2b[good_evt]->SetTitle(Form("dtheta 12C vs 16O vs 16O %g %g %g",dtheta[0],dtheta[1],dtheta[2]));
			if(good_evt<draw_events)hdistxy2c[good_evt]->SetTitle(Form("dL 12C light heavy %g %g",dlength[0][0],dlength[0][1]));
			if(good_evt<draw_events)hdistyz2c[good_evt]->SetTitle(Form("dL 16O light heavy %g %g",dlength[1][0],dlength[1][1]));
			if(good_evt<draw_events)hdistxz2c[good_evt]->SetTitle(Form("dL 16O* light heavy %g %g",dlength[2][0],dlength[2][1]));
			good=true;
			if(fabs(fit_theta[0][0]-fit_theta[0][1])<1 && fabs(fit_theta[1][0]-fit_theta[1][1])<1) {
				good=false;
				cout<<"Straight line"<<endl;
				continue;
			}
//			good=false;
			good=true;
			treetrack->Fill();

			double smallest_chi=min(chisqrd16O,chisqrd16Op);
			smallest_chi=min(smallest_chi,chisqrd12C);
			if(smallest_chi>chisqrdcut) continue;
			if(chisqrd12C<chisqrd16O && chisqrd12C<chisqrd16Op) hchi12Ctheta->Fill(chisqrd12C,fit_theta[0][0]);
			if(chisqrd16O<chisqrd12C && chisqrd16O<chisqrd16Op) hchi16Otheta->Fill(chisqrd16O,fit_theta[1][0]);
			if(chisqrd16Op<chisqrd16O && chisqrd16Op<chisqrd12C) hchi16Optheta->Fill(chisqrd16Op,fit_theta[2][0]);

			if(chisqrd16O>chisqrdcut && chisqrd12C>chisqrdcut && chisqrd16Op>chisqrdcut) {
				if(good_evt<draw_events) {
					htrackxz[good_evt]->SetTitle("Bad chisqrd");
					hdistxy2[good_evt]->Reset();
					hdistxy2a[good_evt]->Reset();
					hdistxy2b[good_evt]->Reset();
					hdistxy2c[good_evt]->Reset();
					hdistyz2[good_evt]->Reset();
					hdistyz2a[good_evt]->Reset();
					hdistyz2b[good_evt]->Reset();
					hdistyz2c[good_evt]->Reset();
					hdistxz2[good_evt]->Reset();
					hdistxz2a[good_evt]->Reset();
					hdistxz2b[good_evt]->Reset();
					hdistxz2c[good_evt]->Reset();
				}
				Ransac2();//Lets try again
				select_nalpha=1;
				double take2_12C=Fitnalpha();
				chisqrd12C=take2_12C;
				select_nalpha=2;
				double take2_16O=Fitnalpha();
				chisqrd16O=take2_16O;
				select_nalpha=3;
				double take2_16Op=Fitnalpha();
				chisqrd16Op=take2_16Op;

				if(chisqrd16O>chisqrdcut && chisqrd12C>chisqrdcut && chisqrd16Op>chisqrdcut) {
					if(good_evt<draw_events) htrackxz[good_evt]->SetTitle("Bad chisqrd - 2 tries");
					good=true;
					continue;
				}
			}
			if(length[0]>length[1]) hthetathetaraw->Fill(fit_theta[0][0],fit_theta[0][1]);
			if(length[0]<=length[1]) hthetathetaraw->Fill(fit_theta[0][1],fit_theta[0][0]);
			if(good_evt<draw_events)hdistxz2a[good_evt]->SetTitle(Form("chisqrd (12C/16O/16O*) %g %g %g",chisqrd12C,chisqrd16O,chisqrd16Op));

			if(vertexy<20 || vertexy>160) {
				if(good_evt<draw_events) htrackxz[good_evt]->SetTitle("Vertex removed");
				good=true;
				continue;//Remove vertices too close to the edge!!
			}
			hmultE2->Fill(mult_decay,peakEnergy);
			hselect->Fill((chisqrd12C-chisqrd16O)/(chisqrd12C+chisqrd16O));
			hselect12C16O->Fill((chisqrd12C-chisqrd16O)/(chisqrd12C+chisqrd16O));
			hselect12C16Op->Fill((chisqrd12C-chisqrd16Op)/(chisqrd12C+chisqrd16Op));
			hselect16O16Op->Fill((chisqrd16O-chisqrd16Op)/(chisqrd16O+chisqrd16Op));
//			hselect->Fill((chisqrd12C-chisqrd16O)/(chisqrd12C+chisqrd16O));
			hselect3->Fill((chisqrd12C-chisqrd16O)/(chisqrd12C+chisqrd16O),vertexy);
			hmultspread2->Fill(mult_decay,spread);
			double smallest_chi2=min(chisqrd216O,chisqrd216Op);
			hchisqrdnalpha->Fill(smallest_chi);
			hconf212C->Fill((chisqrd212C-smallest_chi2)/(chisqrd212C+smallest_chi2));
			smallest_chi2=min(chisqrd216Op,chisqrd212C);
			hconf216O->Fill((chisqrd216O-smallest_chi2)/(chisqrd216O+smallest_chi2));
			smallest_chi2=min(chisqrd216O,chisqrd212C);
			hconf216Op->Fill((chisqrd216Op-smallest_chi2)/(chisqrd216Op+smallest_chi2));

			if(good_evt<draw_events) htrackxz[good_evt]->SetTitle("Uncategorized");
//			if(dtheta[0]<=dtheta[1]){
			double min_chi=min(chisqrd16O,chisqrd16Op);
			hconf12C->Fill((chisqrd12C-min_chi)/(chisqrd12C+min_chi));
			hconf12C2D->Fill((chisqrd12C-min_chi)/(chisqrd12C+min_chi),chisqrd12C);
			min_chi=min(chisqrd12C,chisqrd16Op);
			hconf16O->Fill((chisqrd16O-min_chi)/(chisqrd16O+min_chi));
			hconf16O2D->Fill((chisqrd16O-min_chi)/(chisqrd16O+min_chi),chisqrd16O);
			min_chi=min(chisqrd12C,chisqrd16O);
			hconf16Op->Fill((chisqrd16Op-min_chi)/(chisqrd16Op+min_chi));
			hconf16Op2D->Fill((chisqrd16Op-min_chi)/(chisqrd16Op+min_chi),chisqrd16Op);
			if(chisqrd12C<chisqrd16O && chisqrd12C<chisqrd16Op) {
				hthetatheta->Fill(fit_theta[0][0],fit_theta[0][1]);
				scattering_theta[0]=fit_theta[0][0];
				scattering_theta[1]=fit_theta[0][1];
				hheavyR->Fill(rheavy[0],rtheavy[0]);
				hlightR->Fill(rlight[0],rtlight[0]);
				good=true;
				h3_g0->Fill(totalEnergy);

//				if(dtheta[0]<10.) {
					hthetatheta12C->Fill(fit_theta[0][0],fit_theta[0][1]);
					htheta12C->Fill(fit_theta[0][0]);
//					cout<<"lab/cm: "<<fit_theta[0][0]<<"\t"<<GetCMAngle12C(fit_theta[0][0])<<endl;
					htheta12CCM->Fill(GetCMAngle12C(fit_theta[0][0]));
					if(good_evt<draw_events) htrackxz[good_evt]->SetTitle(Form("12C(n,a0) %g %g",fit_theta[0][0],fit_theta[0][1]));

					cout<<"12C(n,a0)"<<endl;
//				}
			}
//			if(dtheta[0]>dtheta[1]){
			if(chisqrd16O<=chisqrd12C && chisqrd16O<chisqrd16Op) {
				hthetatheta->Fill(fit_theta[1][0],fit_theta[1][1]);
				scattering_theta[0]=fit_theta[1][0];
				scattering_theta[1]=fit_theta[1][1];
				hheavyR->Fill(rheavy[1],rtheavy[1]);
				hlightR->Fill(rlight[1],rtlight[1]);
				good=true;
				h3_g1->Fill(totalEnergy);

//				if(dtheta[1]<10.) {
					hthetatheta16O->Fill(fit_theta[1][0],fit_theta[1][1]);
					htheta16O->Fill(fit_theta[1][0]);
					htheta16OCM->Fill(GetCMAngle16O(fit_theta[1][0]));
//					cout<<"16O(n,a0)"<<endl;
//					cout<<"lab/cm: "<<fit_theta[1][0]<<"\t"<<GetCMAngle12C(fit_theta[1][0])<<endl;
					double min_chi=min(chisqrd12C,chisqrd16Op);

					if(good_evt<draw_events) htrackxz[good_evt]->SetTitle(Form("16O(n,a0) %g %g",fit_theta[1][0],fit_theta[1][1]));
//				}
			}
			if(chisqrd16Op<=chisqrd12C && chisqrd16Op<chisqrd16O) {
				hthetatheta->Fill(fit_theta[2][0],fit_theta[2][1]);
				scattering_theta[0]=fit_theta[2][0];
				scattering_theta[1]=fit_theta[2][1];
				hheavyR->Fill(rheavy[2],rtheavy[2]);
				hlightR->Fill(rlight[2],rtlight[2]);
				h3_g2->Fill(totalEnergy);

//				if(dtheta[2]<10.) {
					hthetatheta16Op->Fill(fit_theta[2][0],fit_theta[2][1]);
					htheta16Op->Fill(fit_theta[2][0]);
//					cout<<"lab/cm: "<<fit_theta[2][0]<<"\t"<<GetCMAngle12C(fit_theta[2][0])<<endl;

					htheta16OpCM->Fill(GetCMAngle16O(fit_theta[2][0]));
//					cout<<"16O(n,a1)"<<endl;
					double min_chi=min(chisqrd16O,chisqrd12C);
					if(good_evt<draw_events) htrackxz[good_evt]->SetTitle(Form("16O(n,a1) %g %g",fit_theta[2][0],fit_theta[2][1]));


//				}
			}
			h3_g->Fill(totalEnergy);

//			hthetatheta->Fill(max(fit_theta[0][0],fit_theta[0][1]),min(fit_theta[0][0],fit_theta[0][1]));
//			Fitnalpha16O();
//			HoughTransform(1);//find a
//			HoughTransform(2);//find 13C
			if(length[0]>length[1]) {//assumes alpha is longer == length[0]
				neutronE16O=GetEnergy4He(length[0])+GetEnergy13C(length[1]);
				energy[0]=GetEnergy4He(length[0]);
				energy[1]=GetEnergy13C(length[1]);
				neutronE12C=GetEnergy4He(length[0])+GetEnergy13C(length[1]);
			}
			if(length[1]>=length[0]) {
				neutronE16O=GetEnergy4He(length[1])+GetEnergy13C(length[0]);
				energy[0]=GetEnergy4He(length[1]);
				energy[1]=GetEnergy13C(length[0]);
				neutronE12C=GetEnergy4He(length[1])+GetEnergy9Be(length[0]);
			}
			good=true;//
//			cout<<"Assuming n,a "<<length[0]<<"\t"<<length[1]<<endl;
			if(good_evt<draw_events)hdistxy2a[good_evt]->SetTitle(Form("Thetas (%g %g) (%g %g) (%g %g)",fit_theta[0][0],fit_theta[0][1],fit_theta[1][0],fit_theta[1][1],fit_theta[2][0],fit_theta[2][1]));
//			if(good_evt<draw_events)hdistyz2a[good_evt]->SetTitle(Form("Ranges  %g %g, Theta %g %g",length[0],length[1],scattering_theta[0],scattering_theta[1]));
//			if(good_evt<draw_events)hdistxz2a[good_evt]->SetTitle(Form("Qvalue (12C/16O) %g %g",neutronE12C,neutronE16O));
/*
			if((dtheta[0]<5 || dtheta[1]<5) && (chisqrd12C<10 || chisqrd16O<10)) {//
//Q-value = -5.702052 +/- 0.000077 MeV        for the reaction: 12C(1n,4He)9Be
//Q-value = -2.215609 +/- 0.000001 MeV        for the reaction: 16O(1n,4He)13C
				cout<<"Reasonable n,a "<<evtno<<endl;
				if(fabs(scattering_theta[0]+scattering_theta[1]-180)<10) {
//					good=false;
//					continue;//this is just a straight line that has been fitted as two halves
				}
				hQvalue12C->Fill(neutronE12C);//Total energy assuming this was 12C(n,a) from lengths
				hQvalue16O->Fill(neutronE16O);//Total energy assuming this was 16O(n,a) from lengths
				hQvalue2D->Fill(neutronE12C,neutronE16O);//
//				cout<<"n,a  "<<neutronE12C<<"\t"<<neutronE16O<<"\t"<<length[0]<<"\t"<<length[1]<<"\t"<<GetEnergy4He(length[1])<<"\t"<<GetEnergy9Be(length[0])<<endl;
//				hthetatheta->Fill(scattering_theta[0],scattering_theta[1]);
				hEavg->Fill(arm_energy[0]/length[0],arm_energy[1]/length[1]);
//				if(good_evt<draw_events)hdistyz2a[good_evt]->SetTitle(Form("Eavg) %g %g",arm_energy[0]/length[0],arm_energy[1]/length[1]));
				good=false;
				good=true;
				if(good_evt<draw_events) htrackxy[good_evt]->SetTitle(Form("n,a theta theta: %g %g",scattering_theta[0],scattering_theta[1]));

			}
			else {
				fexperiment=1;//If shorter then must be a Hoyle event
				good=true;
			}
*/
		}
		if(fexperiment==1) {//n,n2
			good=true;
			if(verbose)cout<<"Analysing Hoyle event"<<endl;
			double Ecarbon=0,Ealphas=0;
//			for(int k=1;k<4;k++) HoughTransform(k);//find 3 lines
			int goodlines[3]={0};
//			for(int k=0;k<3;k++) goodlines[k]=HoughTransformEvent();
			Ransac();
			Ransac3();
			good=true;
			if(verbose)cout<<"YZ "<<goodlines[0]<<"\t"<<goodlines[1]<<"\t"<<goodlines[2]<<endl;

			double pcarbon[4]={0.};
			for(int k=0;k<3;k++) {//loop over alphas
				energy[k]=GetEnergy(length[k]);
//				energy[k]=GetEnergy4He(length[k]);
//				cout<<"ENERGY: "<<k<<"\t"<<energy[k]<<endl;
				for(int m=0;m<3;m++) {//loop over xyz
					if(energy[k]>0) {
						palpha[k][m]/=length[k];//normalize the lengths first
						palpha[k][m]*=sqrt(2.*4.*energy[k]);//times by p_total
						pcarbon[m]+=palpha[k][m];//sum of all alphas is carbon momentum in xyz
					}
					else {
						palpha[k][m]=0.;
					}
				}
			}
			hlinearity->Fill(chi_linear,chi_triple);
//			hChi1->Fill(chi_linear);
//			hChi2->Fill(chi_triple);
			bool find3lines=true;
			if((chi_linear>6. && chi_triple<4) || chi_linear>10.) find3lines=true;
//			cout<<"CHI linear/triple: "<<chi_linear<<"\t"<<chi_triple<<" for event "<<evtno<<endl;
///			if(goodlines[0]==1 && goodlines[1]==1 && goodlines[2]==1) find3lines=true;
//			if((energy[0]==0 || energy[1]==0 || energy[2]==0)) find3lines=false;
//			cout<<"energies: "<<energy[0]<<"\t"<<energy[1]<<"\t"<<energy[2]<<endl;
//			cout<<"DELIBERATELY SELECTING BAD EVENTS"<<endl;
			if(!find3lines) {
				fexperiment=2;//Probably just an n,a
				if(verbose)cout<<"L: "<<length[0]<<"\t"<<length[1]<<"\t"<<length[2]<<endl;
				if(verbose)cout<<"E: "<<energy[0]<<"\t"<<energy[1]<<"\t"<<energy[2]<<endl;
				if(verbose)cout<<"Just n,n0 "<<goodlines[0]<<"\t"<<goodlines[1]<<"\t"<<goodlines[2]<<endl;
				good=false;
//				good=true;
//				fexperiment=1;
//				continue;//
			}
			else {
				good=true;
				Ecarbon=0.;
				Ealphas=0.;
				for(int m=0;m<3;m++) {//xyz
					Ecarbon+=pcarbon[m]*pcarbon[m];//p^2
					Ealphas+=energy[m];
				}
				pcarbon[3]=sqrt(Ecarbon);//Ecarbon here is p^2, pcarbon[3] is p_total
				Ecarbon/=(2.*12);//E=p^2/(2*m)
				if(verbose)for(int k=0;k<3;k++) cout<<"Alpha "<<k<<"\t"<<palpha[k][0]<<"\t"<<palpha[k][1]<<"\t"<<palpha[k][2]<<"\t"<<energy[k]<<"\t"<<length[k]<<endl;
				for(int k=0;k<3;k++) cout<<"Alpha p[xyz],E,l "<<k<<"\t"<<palpha[k][0]<<"\t"<<palpha[k][1]<<"\t"<<palpha[k][2]<<"\t"<<energy[k]<<"\t"<<length[k]<<endl;
				if(good_evt<draw_events) htrackxy[good_evt]->SetTitle(Form("alpha energies for Hoyle: %g %g %g",energy[0],energy[1],energy[2]));
				if(good_evt<draw_events) htrackyz[good_evt]->SetTitle(Form("Ex Hoyle: %g",Ealphas-Ecarbon+7.272));
				if(good_evt<draw_events) htrackxz[good_evt]->SetTitle(Form("Length alphas:%g %g %g",length[0],length[1],length[2]));
				good=true;
				hEx->Fill(Ealphas-Ecarbon+7.272);//Ex = sum_i E_alpha_i - E_carbon - Q
			cout<<"Ex: "<<Ealphas-Ecarbon+7.272<<"\t"<<Ecarbon<<endl;
				hErEs->Fill(Ealphas,totalEnergy);
//			cout<<totalEnergy<<endl;
				hEx_d->Fill(Ealphas-Ecarbon+7.272,maxy);
				Ex=Ealphas-Ecarbon+7.272;
				if(verbose)cout<<"3 alphas: "<<energy[0]<<"\t"<<energy[1]<<"\t"<<energy[2]<<endl;
				if(good_evt<draw_events) hdistxz3[good_evt]->SetTitle(Form("Range E, WF E: %g %g",energy[0]+energy[1]+energy[2],totalEnergy));
//			cout<<"SCALE: "<<Ecarbon<<"\t"<<totalEnergy<<endl;
//				Ecarbon=(totalEnergy-45000)/122049;
				double thetaC12=acos(pcarbon[1]/pcarbon[3]);//theta=acos(z/r) --> theta=acos(y/r) in TexAT coordinates
//			hEtheta->Fill(thetaC12/d2r,Ecarbon);
				double palphacm[3][3],Erel[3],dalitzEtotal=0;
//
//			p_n=(p_c^2 + p_c^2/12 - 2*Ex)/(2*p_c*cos(theta))
//
				double p_n,E_n,p_c;
				p_c=sqrt(2.*12*Ecarbon);
				p_n=(p_c*p_c*(13./12.)+2*7.654)/(2.*p_c*cos(thetaC12));
				E_n=(p_n*p_n*0.5);
//			cout<<p_c<<"\t"<<Ecarbon<<"\t"<<thetaC12/d2r<<"\t"<<E_n<<endl;
				hneutronenergy->Fill(E_n);
				hneutronenergy2D->Fill(E_n,thetaC12/d2r);
				if(good_evt<draw_events) hdistxy3[good_evt]->SetTitle(Form("Neutron E, angle = %g %g",E_n,thetaC12/d2r));
				for(int k=0;k<3;k++) {//alphas
					for(int m=0;m<3;m++) {//xyz
						palphacm[k][m]=palpha[k][m]-pcarbon[m];
//					cout<<k<<"\t"<<m<<"\t"<<palphacm[k][m]<<"\t"<<palpha[k][m]<<"\t"<<pcarbon[m]<<endl;
					}
					Erel[k]=(palphacm[k][0]*palphacm[k][0]+palphacm[k][1]*palphacm[k][1]+palphacm[k][2]*palphacm[k][2])/(2.*4.);
					dalitzEtotal+=Erel[k];
				}
				int largest=-1,choice1=-1,choice2=-1;
//			cout<<Erel[0]<<"\t"<<Erel[1]<<"\t"<<Erel[2]<<endl;
				if(Erel[2]>Erel[1] && Erel[2]>Erel[0]) {
					largest=2;
					choice1=0;
					choice2=1;
				}
				if(Erel[1]>Erel[2] && Erel[1]>Erel[0]) {
					largest=1;
					choice1=0;
					choice2=2;
				}
				if(Erel[0]>Erel[1] && Erel[0]>Erel[2]) {
					largest=0;
					choice1=1;
					choice2=2;
				}
//			Reconstruct 8Be
				double PBe[4]={0.},PBet=0.;
				for(int k=0;k<3;k++) {
					PBe[k]=palpha[choice1][k]+palpha[choice2][k];
					PBet+=PBe[k];
				}
				PBe[3]=sqrt(PBet);//total momentum
				double EBe=PBet/(2.*8.);
				hBe->Fill(energy[choice1]+energy[choice2]-EBe);
//			cout<<"CHOICE "<<choice1<<"\t"<<choice2<<endl;
//			cout<<"ENERGIES OF CHOICES "<<energy[choice1]<<"\t"<<energy[choice2]<<"\t"<<EBe<<endl;
				for(int k=0;k<3;k++) Erel[k]/=dalitzEtotal;
				sort(Erel,Erel+3);
				hDalitz->Fill((1./sqrt(3))*(Erel[1]-Erel[0]),(1./3.)*(Erel[1]+Erel[0]-2*Erel[2]));
				}
			}

		if(fexperiment==3) {//n,n0 scattering
//			HoughTransform(1);//find 1 line
			good=false;
			Ransac();
			double carbon[4];
//                        double C12energy=GetEnergy12C(length);
			for(int m=0;m<3;m++) carbon[m]=palpha[0][m];
			double thetaC12,Ecarbon;
			carbon[3]=sqrt(carbon[0]*carbon[0]+carbon[1]*carbon[1]+carbon[2]*carbon[2]);
			thetaC12=acos(carbon[1]/carbon[3]);
//			Ecarbon=totalEnergy*0.000005;
			hEtheta->Fill(scatterangle/d2r,length[0]);
			hscatterEl->Fill(peakEnergy,length[0]);
			if(good_evt<draw_events)hdistxy1[good_evt]->SetTitle(Form("angle: %g   length: %g",scatterangle/d2r,length[0]));
			if(good_evt<draw_events) htrackxy[good_evt]->SetTitle(Form("angle: %g   length: %g",scatterangle/d2r,length[0]));
//			cout<<"Assuming n,n0"<<endl;
		}

/*

		if(fexperiment==1) {
			good=true;
			if(verbose)cout<<"Interesting event? "<<good<<endl;
		}
		else {
			good=false;
		}
*/
		if(debug) cout<<"Hough done"<<endl;
		if(readmc) {
			beamx=66;
			beamy=56;
			beamz=50;
		}
		hcount->Fill(fexperiment);
                hcount->GetXaxis()->SetBinLabel(1,"Hoyle");
                hcount->GetXaxis()->SetBinLabel(2,"n,a");
                hcount->GetXaxis()->SetBinLabel(3,"n,n0");

                hevents->Fill(event_selection);
                hevents->GetXaxis()->SetBinLabel(event_selection+1,"FOUND A GOOD CENTER");
                event_selection++;
///////////////////////////////////////////////////////////////////////////////////////
////		Check distance between beam end and the found center of the decay ////
//////////////////////////////////////////////////////////////////////////////////////
		if(debug) cout<<"Fill rndm walk"<<endl;
		double centr_dist_x,centr_dist_y,centr_dist_z,beam_dist_x,beam_dist_y,beam_dist_z;
//		Get positions
		Dist(maxx,maxy,maxz,centr_dist_x,centr_dist_y,centr_dist_z);
////////////////////////////////////////////////////////////
//		Find furthest point
////////////////////////////////////////////////////////////
		bool good_arm=true;
		bool doubleup=false;
		double arm[3][2]={0.};
		int npeaks_ang_=-1;
		bool forward_arm=false;
		int tries=0,tries_center=0;
		bool moved_center=false;
                hevents->Fill(event_selection);
                hevents->GetXaxis()->SetBinLabel(event_selection+1,"ANOTHER ENERGY CUT");
                event_selection++;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(good_evt<draw_events) {
//			htrackxy[good_evt]->Fill(dmaxx,dmaxy,1000000);//mark these points with a large weighted bin
			double x_start,y_start,z_start;
			Dist(dmaxx,dmaxy,dmaxz,x_start,y_start,z_start);
//			hdistxy[good_evt]->Fill(x_start,y_start,1000000);//mark these points with a large weighted bin
//			hdistxz[good_evt]->Fill(x_start,z_start,1000000);//mark these points with a large weighted bin
//			hdistyz[good_evt]->Fill(y_start,z_start,1000000);//mark these points with a large weighted bin
//			htrackyz[good_evt]->Fill(dmaxy,dmaxz,1000000);
//			htrackxy[good_evt]->Fill(dmaxx,dmaxy,2000000);
//			htrackyz[good_evt]->Fill(dmaxy,dmaxz,2000000);
//			htrackxz[good_evt]->Fill(dmaxx,dmaxz,2000000);
//			htrackxz[good_evt]->Fill(distx,distz,2000000);
//			if(fexperiment==2)htrackxy[good_evt]->SetTitle(Form("Track lengths=%g %g %g",length[0],length[1],length[2]));//record the total energy on the 3D track histogram
//			if(fexperiment==2)htrackyz[good_evt]->SetTitle(Form("Total energy=%g",neutronE));
			if(fexperiment==1)hdistxy[good_evt]->SetTitle(Form("Track lengths=%g %g %g",length[0],length[1],length[2]));//record the total energy on the 3D track histogram
			hdistyz[good_evt]->SetTitle(Form("Ex=%g,Spreadiness=%g",Ex,largest_width));
			if(fexperiment==1)hdistxz[good_evt]->SetTitle(Form("Energies=%g %g %g",energy[0],energy[1],energy[2]));
			//hdistxy1[good_evt]->SetTitle(Form("Chi_linear vs chi_3 %g %g",chi_linear,chi_triple));
//			if(fexperiment==2)htrackxz[good_evt]->SetTitle(Form("Energies=%g %g",energy[0],energy[1]));
/*
			for(int l=0;l<3;l++) {
				if(goodarms){
					htrackxy[good_evt]->Fill(armxyz[l][0],armxyz[l][1],1500000);
					htrackyz[good_evt]->Fill(armxyz[l][1],armxyz[l][2],1500000);
					htrackxz[good_evt]->Fill(armxyz[l][0],armxyz[l][2],1500000);
					if(debug)cout<<"ARMXYZ\t"<<armxyz[l][0]<<"\t"<<armxyz[l][1]<<"\t"<<armxyz[l][2]<<endl;
				}
			}
*/
//			htrack3D[good_evt]->SetTitle(Form("Energy=%g",(peakEnergy)));//record the total energy on the 3D track histogram
		}
		if(good || takeall || readmc || !twopmode) {
			hEtime->Fill(peakEnergy,mmD2PTime);
			for(int m=0;m<6;m++) {
				if(firstpadE[m]>0)Ebeam+=firstpadE[m];//find total value
			}
//			hdeE->Fill(Ebeam,peakEnergy);
			hdeL->Fill(Ebeam,beamy);
			if(!mixedevent)hbeamtime->Fill(Ebeam,mmD2PTime);//check to see if the energy left from the beam decays over time in the memory array
			hbeamEdE->Fill(Ebeam,beamy);
		}
	}
	good=true;
	delete tree;
	return;
}


void Analysis::Deconvolution() {
//		cout<<"DECONV "<<debug<<endl;
		TSpectrum *s2 = new TSpectrum();
//		for(int t=0;t<500;t++) cout<<"DECONV: "<<t<<"\t"<<wave[t]<<"\t"<<response[t]<<endl;
		const char * errorval = s2->Deconvolution(wave,response,timetemp,15,15,1);//deconvolute
		if(errorval!=0) cout<<"Deconv Error!: "<<errorval<<endl;
		for(int t=0;t<timetemp;t++) {
//			if(t<200) wave[t]=0;
//			if(wave[t]<0 || t<source_offset) wave[t]=0.;
			hwave->SetBinContent(t,wave[t]);
//			cout<<"RESULT: "<<t<<"\t"<<wave[t]<<endl;
//			else hwave->SetBinContent(t,0);
			if(wave[t]<200) {
				hwave->SetBinContent(t,0.);
				wave[t]=0;
			}
			if(good_evt<draw_events) hwaveform_dc[good_evt]->Fill(waveformno,t,wave[t]);
		}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////// /// Get the multiple? peaks from the deconvoluted spectrum
		const int npeaks = s2->Search(hwave,5,"",0.1);//peak finder, sigma val, threshold of 0.1 (if peak is 0.1*maxvalue then count)
		double *xpeaks = s2->GetPositionX();//list of peaks
		const int good_peaks=good_double;
		good_double=0;//reset
		if(debug) cout<<"Waveform: "<<waveformno<<" found "<<npeaks<<" peaks"<<endl;
//		cout<<"Waveform: "<<waveformno<<" found "<<npeaks<<" peaks"<<endl;
		for(int l=0;l<npeaks;l++) {//If peak is above 300 then take it as a "Proper" peak
			if(debug) cout<<xpeaks[l]<<"\t"<<hwave->GetBinContent(xpeaks[l])<<endl;
//			if((hwave->GetBinContent(xpeaks[l])>180. && xpeaks[l]>source_offset+30 && !readmc && xpeaks[l]<timetemp-200) || (mmWaveformY[waveformno][int(xpeaks[l]-source_offset)]>mc_thresh && readmc && xpeaks[l]>200 && xpeaks[l]<450)) {
//			cout<<l<<"\t of "<<npeaks<<"\t"<<xpeaks[l]<<"\t"<<hwave->GetBinContent(xpeaks[l]+1)<<endl;
			if((hwave->GetBinContent(xpeaks[l]+1)>300 && xpeaks[l]>20 && xpeaks[l]<timetemp-100 && (hwave->GetBinContent(xpeaks[l])>50 || hwave->GetBinContent(xgoodpeaks[l]+2)>50))) {
				xgoodpeaks[good_double]=xpeaks[l];
				if(debug)cout<<"Peak "<<l<<"->"<<good_double<<" is "<<xpeaks[l]<<endl;
//				cout<<"Peak "<<l<<"->"<<good_double<<" is "<<xpeaks[l]<<"\t"<<hwave->GetBinContent(xpeaks[l])<<endl;
				good_double++;
			}
		}
		if(good_double==0 && !readmc) {//Something has gone awry with the deconvolution and extraction
			double maxE=0.,maxT=0.;
			for(int t=0;t<timetemp;t++) { // if(debug)cout<<"TEST\t"<<t<<"\t"<<mmWaveformY[j][t]<<endl;
				if(mmWaveformY[waveformno][t]>maxE && t>source_offset && t<400) {
					maxE=mmWaveformY[waveformno][(t-source_offset+512)%512];
					maxT=t;
				}
		               	if(maxE>150. || (readmc && maxE>mc_thresh)) {
					xgoodpeaks[0]=maxT-source_offset;
					good_double=1;
				}
			}
			if(debug)cout<<"SINGLE: "<<xgoodpeaks[0]<<endl;
		}
		sort(xgoodpeaks,xgoodpeaks+good_double);//order
//		cout<<good_double<<" good hits"<<endl;

return;
}

void Analysis::WriteHist(TString histname) { //	Add histograms to TList which we can then write to a ROOT file. Options to turn off certain options are defined below (if crashing due to large file size)
	cout<<"Writing histograms\t"<<histname<<endl;
	TList *list = new TList();
	list->Add(hevents);
	list->Add(hcount);
	list->Add(h1);
	list->Add(h2);
	list->Add(hErEs);
	list->Add(hE1E2);
	list->Add(helscat);
	list->Add(hlinearity);
	list->Add(hChi1);
	list->Add(hChi2);
	list->Add(h3);
	list->Add(h3_g);
	list->Add(h3_g0);
	list->Add(h3_g1);
	list->Add(h3_g2);
	list->Add(hdEE);
	list->Add(hdeL);
	list->Add(hEx_d);
	list->Add(hEx_evt);
	list->Add(hEx1Ex2);
	list->Add(hpeaktotal);
	list->Add(hbeamenergy);
	list->Add(hEtime);
	list->Add(hheavyR);
	list->Add(hlightR);
	hthetatheta->SetTitle(";Light angle (degrees);Heavy angle (degrees)");
	list->Add(hthetatheta);
	list->Add(hthetathetaraw);
	hthetathetatheory->SetTitle(";Light angle (degrees);Heavy angle (degrees)");
	list->Add(hthetatheta12C);
	list->Add(hthetatheta16O);
	list->Add(hthetatheta16Op);
	list->Add(hthetathetatheory);
	list->Add(hthetathetatheoryCM);
	htheta12C->SetTitle(";Light angle (degrees);Counts");
	list->Add(hselect);
	list->Add(hselect12C16O);
	list->Add(hselect12C16Op);
	list->Add(hselect16O16Op);
	list->Add(hconf12C);
	list->Add(hconf16O);
	list->Add(hconf16Op);
	list->Add(hconf212C);
	list->Add(hconf216O);
	list->Add(hconf216Op);
	list->Add(hconf12C2D);
	list->Add(hconf16O2D);
	list->Add(hconf16Op2D);
	list->Add(hchi12Ctheta);
	list->Add(hchi16Otheta);
	list->Add(hchi16Optheta);
	list->Add(hchisqrdth12C);
	list->Add(hchisqrdth16O);
	list->Add(hchisqrdth16Op);
	list->Add(hchisqrdE12C);
	list->Add(hchisqrdE16O);
	list->Add(hchisqrdE16Op);
	list->Add(hselect2);
	list->Add(hselect3);
	list->Add(htheta12C);
	htheta16O->SetTitle(";Light angle (degrees);Counts");
	list->Add(htheta16O);
	list->Add(htheta16Op);
	htheta12CCM->SetTitle(";Light angle cm (degrees);Counts");
	list->Add(htheta12CCM);
	htheta16OCM->SetTitle(";Light angle cm(degrees);Counts");
	list->Add(htheta16OCM);
	list->Add(htheta16OpCM);
	list->Add(hchisqrd);
	list->Add(hchisqrdnalpha);
	list->Add(hEavg);
	list->Add(hwidth);
	list->Add(htime);
	list->Add(htime2);
	list->Add(hallhits);
	list->Add(hallhits2);
	list->Add(hEx);
	list->Add(hneutronenergy);
	list->Add(hneutronenergy2D);
	list->Add(hpeaks);
	list->Add(hDalitz);
	list->Add(hEtheta);
	list->Add(hEthetatheory);
	list->Add(hBe);
	list->Add(hmult);
	list->Add(hmultspread);
	list->Add(hmultspread2);
	list->Add(hmultE);
	list->Add(hmultE2);
	list->Add(hscatterEl);
	list->Add(hQvalue12C);
	list->Add(hQvalue16O);
	hQvalue2D->SetTitle("12C vs 16O q-value (-5.7 vs -2.2)");
	list->Add(hQvalue2D);
	for(int i=0;i<20;i++) {
//		list->Add(hSilicon[i]);
//		list->Add(hSilicong[i]);
//		list->Add(hSilicongg[i]);
	}
	bool write_by_event=true;
	if(write_by_event && !server_mode) {
//		cout<<good_evt<<"\tgood events - reading"<<min(draw_events,good_evt)<<endl;
		for(int i=0;i<min(draw_events,good_evt);i++) {
//			cout<<i<<"\t"<<draw_events<<endl;
//			list->Add(hpulse[i]);
//			list->Add(hpulse_b[i]);
//			list->Add(hwaveform[i]);
//			list->Add(hwaveform_dc[i]);
//			list->Add(hwaveform_hit[i]);
//			list->Add(hHits[i]);
//			list->Add(hHits_side[i]);
//			list->Add(hpulse_chn[i]);
//			list->Add(hpulse_chndc[i]);
//			list->Add(htrack3D[i]);
//			list->Add(henergy[i]);
			list->Add(htrackxy[i]);
			list->Add(htrackyz[i]);
			list->Add(htrackxz[i]);
//			list->Add(hflat[i]);
//			list->Add(hbeam[i]);
//			list->Add(hbeam2[i]);
//			list->Add(htrackxysmooth[i]);
//			list->Add(hdeltat[i]);
//			list->Add(hdot[i]);
//			list->Add(hdedx[i]);
//			list->Add(hangle_arms[i]);
//			list->Add(hsuper[i]);
//			list->Add(hnearestneighbour[i]);
//			list->Add(harm1_xy[i]);
//			list->Add(harm1_yz[i]);
//			list->Add(hmults[i]);
//			list->Add(hHoughXY[i]);
//			list->Add(hHoughYZ[i]);
//			list->Add(hHoughXZ[i]);
//			list->Add(hfitbeamXY[i]);
//			list->Add(hfitbeamYZ[i]);
//			list->Add(hHoughXY2[i]);
//                        list->Add(hHoughYZ2[i]);
//                      list->Add(hHoughXZ2[i]);
			list->Add(hdistxy[i]);
			list->Add(hdistyz[i]);
			list->Add(hdistxz[i]);
//			list->Add(hdistxy1[i]);
//			list->Add(hdistyz1[i]);
//			list->Add(hdistxz1[i]);
			list->Add(hdistxy2[i]);
			list->Add(hdistyz2[i]);
			list->Add(hdistxz2[i]);
/*
			list->Add(hdistxy2a[i]);
			list->Add(hdistyz2a[i]);
			list->Add(hdistxz2a[i]);
			list->Add(hdistxy2b[i]);
			list->Add(hdistyz2b[i]);
			list->Add(hdistxz2b[i]);
			list->Add(hdistxy2c[i]);
			list->Add(hdistyz2c[i]);
			list->Add(hdistxz2c[i]);
*/
			list->Add(hdistxy3[i]);
			list->Add(hdistyz3[i]);
			list->Add(hdistxz3[i]);
//			list->Add(harmxy[i]);
//			list->Add(harmyz[i]);
//			list->Add(harmxz[i]);
//			list->Add(hmaxstrips[i]);
//			list->Add(hmaxchains[i]);
//			list->Add(hHoughprojXY[i]);
//			list->Add(hHoughprojYZ[i]);
//			list->Add(hHoughprojXZ[i]);
//			list->Add(hstage1[i]);
//			list->Add(hstage2[i]);
//			list->Add(hstage3[i]);
//			list->Add(hchainsevt[i]);
//			cout<<"Good event "<<i<<" written"<<endl;
		}
	}
	bool writedoublehits=false;
	if(writedoublehits) {
		for(int i=0;i<doublecount;i++) {
		list->Add(hdoublecandidate[i]);
	}
	}
	bool writehist=true;
	if(writehist) {
		TFile *fout = new TFile(histname,"recreate");
		list->Write("myhists",TObject::kSingleKey);
		fout->ls();
		fout->Close();
	}
	return;
}

void Analysis::DrawTheory() {//draw any theory curves
	double En=Ebeam;//MeV neutron energy2
        double xc[1000],yc[1000],xd[1000],yd[1000];
	double mB=1,mC=9,mD=4,Q=-5.70,Ex=0.,theta,psi,x,pi=4.*atan(1.),Ec;
	ofstream outc12na;
	outc12na.open("outc12na.txt");
	outc12na<<"Theta\tpsi\tEBe\tEa"<<endl;
        for(int i=0;i<1000;i++) {//12C(n,a)9Be
                theta=pi*(i/1000.);
                double a,b,c;
                a=(mD/mB+mC/mB);
                b=-2.*sqrt(En*mD/mB)*cos(theta);
                c=(-mC/mB)*(Q-Ex+En)+En;
                x=(-b+sqrt(b*b-4.*a*c))/(2.*a);
                Ec=En+Q-x*x-Ex;
                double Pc=sqrt(2.*mC*Ec);
                double Pd=sqrt(2.*mD*x*x);
//		psi=asin((sqrt(2.*mD*x*x)/sqrt(2.*mC*Ec))*sin(theta));
		psi=asin(sin(theta)*Pd/Pc);
                double Pb=sqrt(2.*mB*En);
                psi=acos((Pb-Pd*cos(theta))/Pc);
                xc[i]=theta/d2r;
                yc[i]=Ec;
                xd[i]=theta/d2r;
                yd[i]=x*x;
		theta12C.push_back(theta/d2r);
		psi12C.push_back(psi/d2r);
		E12C_heavy.push_back(Ec);
		E12C_light.push_back(x*x);
		hthetathetatheory->Fill(theta/d2r,psi/d2r,1);
		outc12na<<theta/d2r<<"\t"<<psi/d2r<<"\t"<<Ec<<"\t"<<x*x<<endl;
		double thetaCM = GetCMAngle12C(theta/d2r);
		hthetathetatheoryCM->Fill(theta/d2r,thetaCM);
//		cout<<theta/d2r<<"\t"<<psi/d2r<<endl;
        }
	outc12na.close();
	Q=-5.7-1.68;
        for(int i=0;i<1000;i++) {//12C(n,a)9Be*(1.68 MeV)
                theta=pi*(i/1000.);
                double a,b,c;
                a=(mD/mB+mC/mB);
                b=-2.*sqrt(En*mD/mB)*cos(theta);
                c=(-mC/mB)*(Q-Ex+En)+En;
                x=(-b+sqrt(b*b-4.*a*c))/(2.*a);
                Ec=En+Q-x*x-Ex;
                double Pc=sqrt(2.*mC*Ec);
                double Pd=sqrt(2.*mD*x*x);
//		psi=asin((sqrt(2.*mD*x*x)/sqrt(2.*mC*Ec))*sin(theta));
		psi=asin(sin(theta)*Pd/Pc);
                double Pb=sqrt(2.*mB*En);
                psi=acos((Pb-Pd*cos(theta))/Pc);
                xc[i]=theta/d2r;
                yc[i]=Ec;
                xd[i]=theta/d2r;
                yd[i]=x*x;
//		theta12C.push_back(theta/d2r);
//		psi12C.push_back(psi/d2r);
//		E12C_heavy.push_back(Ec);
//		E12C_light.push_back(x*x);
//		hthetathetatheory->Fill(theta/d2r,psi/d2r,10);
//		cout<<theta/d2r<<"\t"<<psi/d2r<<endl;
        }
	mC=13.;//16O(n,a)13C
	Q=-2.22;
	ofstream outo16na;
	outo16na.open("outo16na.txt");
	outo16na<<"Theta\tpsi\tEC\tEa"<<endl;

        for(int i=0;i<1000;i++) {
                theta=pi*(i/1000.);
                double a,b,c;
                a=(mD/mB+mC/mB);
                b=-2.*sqrt(En*mD/mB)*cos(theta);
                c=(-mC/mB)*(Q-Ex+En)+En;
                x=(-b+sqrt(b*b-4.*a*c))/(2.*a);
                Ec=En+Q-x*x-Ex;
                double Pc=sqrt(2.*mC*Ec);
                double Pd=sqrt(2.*mD*x*x);
//		psi=asin((sqrt(2.*mD*x*x)/sqrt(2.*mC*Ec))*sin(theta));
		psi=asin(sin(theta)*Pd/Pc);
                double Pb=sqrt(2.*mB*En);
                psi=acos((Pb-Pd*cos(theta))/Pc);
                xc[i]=theta/d2r;
                yc[i]=Ec;
                xd[i]=theta/d2r;
                yd[i]=x*x;

                theta16O.push_back(theta/d2r);
                psi16O.push_back(psi/d2r);
		E16O_heavy.push_back(Ec);
		E16O_light.push_back(x*x);
		hthetathetatheory->Fill(theta/d2r,psi/d2r,100);
		double thetaCM = GetCMAngle16O(theta/d2r);
		hthetathetatheoryCM->Fill(theta/d2r,thetaCM);

		outo16na<<theta/d2r<<"\t"<<psi/d2r<<"\t"<<Ec<<"\t"<<x*x<<endl;
//		cout<<theta/d2r<<"\t"<<psi/d2r<<endl;
        }
	outo16na.close();
	mC=13.;//16O(n,a)13C*(3.09 MeV)
	Q=-2.22-3.09;
        for(int i=0;i<1000;i++) {
                theta=pi*(i/1000.);
                double a,b,c;
                a=(mD/mB+mC/mB);
                b=-2.*sqrt(En*mD/mB)*cos(theta);
                c=(-mC/mB)*(Q-Ex+En)+En;
                x=(-b+sqrt(b*b-4.*a*c))/(2.*a);
                Ec=En+Q-x*x-Ex;
                double Pc=sqrt(2.*mC*Ec);
                double Pd=sqrt(2.*mD*x*x);
//		psi=asin((sqrt(2.*mD*x*x)/sqrt(2.*mC*Ec))*sin(theta));
		psi=asin(sin(theta)*Pd/Pc);
                double Pb=sqrt(2.*mB*En);
                psi=acos((Pb-Pd*cos(theta))/Pc);
                xc[i]=theta/d2r;
                yc[i]=Ec;
                xd[i]=theta/d2r;
                yd[i]=x*x;
                theta16Op.push_back(theta/d2r);
                psi16Op.push_back(psi/d2r);
		E16Op_heavy.push_back(Ec);
		E16Op_light.push_back(x*x);
		hthetathetatheory->Fill(theta/d2r,psi/d2r,1000);
		double thetaCM = GetCMAngle16Op(theta/d2r);
		hthetathetatheoryCM->Fill(theta/d2r,thetaCM);

//		cout<<theta/d2r<<"\t"<<psi/d2r<<endl;
        }
	Q=0;
	mC=12;
	mD=1;
        for(int i=0;i<1000;i++) {//12C(n,el)
                theta=pi*(i/1000.);
                double a,b,c;
                a=(mD/mB+mC/mB);
                b=-2.*sqrt(En*mD/mB)*cos(theta);
                c=(-mC/mB)*(Q-Ex+En)+En;
                x=(-b+sqrt(b*b-4.*a*c))/(2.*a);
                Ec=En+Q-x*x-Ex;
                double Pc=sqrt(2.*mC*Ec);
                double Pd=sqrt(2.*mD*x*x);
//		psi=asin((sqrt(2.*mD*x*x)/sqrt(2.*mC*Ec))*sin(theta));
		psi=asin(sin(theta)*Pd/Pc);
                double Pb=sqrt(2.*mB*En);
                psi=acos((Pb-Pd*cos(theta))/Pc);
		helscat->Fill(-psi/d2r,GetRange12C(Ec));
        }
	Q=0;
	mC=16;
	mD=1;
        for(int i=0;i<1000;i++) {//16O(n,el)
                theta=pi*(i/1000.);
                double a,b,c;
                a=(mD/mB+mC/mB);
                b=-2.*sqrt(En*mD/mB)*cos(theta);
                c=(-mC/mB)*(Q-Ex+En)+En;
                x=(-b+sqrt(b*b-4.*a*c))/(2.*a);
                Ec=En+Q-x*x-Ex;
                double Pc=sqrt(2.*mC*Ec);
                double Pd=sqrt(2.*mD*x*x);
//		psi=asin((sqrt(2.*mD*x*x)/sqrt(2.*mC*Ec))*sin(theta));
		psi=asin(sin(theta)*Pd/Pc);
                double Pb=sqrt(2.*mB*En);
                psi=acos((Pb-Pd*cos(theta))/Pc);
		helscat->Fill(-psi/d2r,GetRange16O(Ec),10);
        }
	Q=-4.44;
	mC=12;
	mD=1;
        for(int i=0;i<1000;i++) {//12C(n,n1)
                theta=pi*(i/1000.);
                double a,b,c;
                a=(mD/mB+mC/mB);
                b=-2.*sqrt(En*mD/mB)*cos(theta);
                c=(-mC/mB)*(Q-Ex+En)+En;
                x=(-b+sqrt(b*b-4.*a*c))/(2.*a);
                Ec=En+Q-x*x-Ex;
                double Pc=sqrt(2.*mC*Ec);
                double Pd=sqrt(2.*mD*x*x);
//		psi=asin((sqrt(2.*mD*x*x)/sqrt(2.*mC*Ec))*sin(theta));
		psi=asin(sin(theta)*Pd/Pc);
                double Pb=sqrt(2.*mB*En);
                psi=acos((Pb-Pd*cos(theta))/Pc);
		helscat->Fill(-psi/d2r,GetRange12C(Ec),100);
        }
	Q=6.05;
	mC=16;
	mD=1;
        for(int i=0;i<1000;i++) {//16O(n,n1)
                theta=pi*(i/1000.);
                double a,b,c;
                a=(mD/mB+mC/mB);
                b=-2.*sqrt(En*mD/mB)*cos(theta);
                c=(-mC/mB)*(Q-Ex+En)+En;
                x=(-b+sqrt(b*b-4.*a*c))/(2.*a);
                Ec=En+Q-x*x-Ex;
                double Pc=sqrt(2.*mC*Ec);
                double Pd=sqrt(2.*mD*x*x);
//		psi=asin((sqrt(2.*mD*x*x)/sqrt(2.*mC*Ec))*sin(theta));
		psi=asin(sin(theta)*Pd/Pc);
                double Pb=sqrt(2.*mB*En);
                psi=acos((Pb-Pd*cos(theta))/Pc);
		helscat->Fill(-psi/d2r,GetRange16O(Ec),1e3);
        }

	return;

}

void Analysis::Dist(int x,int y,int z,double &xout,double &yout,double &zout) {
//		cout<<"DIST IN \t"<<x<<"\t"<<y<<"\t"<<z<<"\t"<<rnd3<<endl;
		double dx=x+rnd3->Rndm();
		double dy=y+rnd3->Rndm();
		double dz=z+rnd3->Rndm();
//		double dx=x;
//		double dy=y;
//		double dz=z;
//		cout<<"RANDS"<<endl;
		if(readmc) dz-=300;
                if(dx<64) xout=dx*1.75+0.5*1.75;
                if(dx>=64 && x<70) xout=63.5*1.75+(dx-63)*3.5;
                if(dx>=70) xout=63.5*1.75+6*3.5+(dx-69)*1.75;
                yout=dy*1.75;
                zout=dz*scaletime;
return;
}


void Analysis::GetCenter(int beamx, int beamy, int beamz, int &maxx, int &maxy, int &maxz,double diff) {
		double maxvalue=0;
		double avg_x,avg_y,avg_z;
		avg_x=hcurrenttrack3D_2->GetMean(1);
		avg_y=hcurrenttrack3D_2->GetMean(2);
		avg_z=hcurrenttrack3D_2->GetMean(3);
		for(int x=64;x<71;x++) {//loop over central pads
                        for(int y=1;y<127;y++) {//loop over length of detector (excluding edge as these have fewer neighbours)
                                for(int z=1;z<timetemp-1;z++) {//loop over the "height" of the detector
					if(hcurrenttrack3D_2->GetBinContent(x+1,y+1,z+1)>0) {
	                                        double distendbeam=sqrt((x-beamx)*(x-beamx)+(y-beamy)*(y-beamy));//ignore Z and offset by 3 to get actual stopping point
        	                                double dist_avg=sqrt((x-avg_x)*(x-avg_x)+(y-avg_y)*(y-avg_y));//distance to mean
                	                        double gauswght=TMath::Exp(-(distendbeam*distendbeam)/(2.*diff*diff));
                        	                double gauswght2=TMath::Exp(-(dist_avg*dist_avg)/(2.*4.*diff*diff));
                                	        if(readmc)gauswght=1.;
                                	        double valuenow=hcurrenttrack3D_2->GetBinContent(x+1,y+1,z+1)*gauswght*gauswght2*hcurrenttrack3D_2->Integral(x,x+2,y,y+2,z,z+2);//To find the "centre" of the decay, use a nearest neighbour interpolation (should hopefully bias to those voxels $
                                	        if(valuenow>maxvalue) {//new max found
                                	                maxx=x;//update
                                	                maxy=y;
                                	                maxz=z;
                                	                maxvalue=valuenow;
                                	        }
//						cout<<x<<"\t"<<y<<"\t"<<z<<"\t"<<hcurrenttrack3D_2->GetBinContent(x,y,z)<<"\t"<<valuenow<<"\t"<<maxvalue<<"\t"<<gauswght2<<"\t"<<hcurrenttrack3D_2->Integral(x-1,x+1,y-1,y+1,z-1,z+1)<<endl;
					}
                                }
                        }
                }
		decay_x=maxx;
		decay_y=maxy;
		decay_z=maxz;
        //maxx,maxy,maxz are the center of the decay
	return;
}

void Analysis::GetEnd(double &maxx, double &maxy, double &maxz) {
	double currentmax=0.;
	double scaling=0.5;
	maxx=-1;
	maxy=-1;
	maxz=-1;
	int ystart=0;
	double x_weight=0,z_weight=0,weight=0;
	bool found=false;
	for(int y=1;y<128;y++) {
		if(found) continue;
		for(int x=63;x<70;x++) {
			for(int z=1;z<timetemp;z++) {
				if(hcurrenttrack3D->GetBinContent(x+1,y+1,z+1)>0.) {//find lowest y with largest value
					maxx=x;
					maxy=y;
					maxz=z;
					x_weight+=1.*x*hcurrenttrack3D->GetBinContent(x+1,y+1,z+1);
					z_weight+=1.*z*hcurrenttrack3D->GetBinContent(x+1,y+1,z+1);
					weight+=1.*hcurrenttrack3D->GetBinContent(x+1,y+1,z+1);
					found=true;
				}
			}
		}
	}
	maxx=(x_weight/weight);
	maxz=(z_weight/weight);
//	cout<<"End: "<<maxx<<"\t"<<maxy<<"\t"<<maxz<<endl;

	return;
}

void Analysis::GetArm(int maxx,int maxy,int maxz,int &distx, int &disty, int &distz) {
	distx=-1;
	disty=-1;
	distz=-1;
	if(maxx==0 || maxy==0 || maxz==0) return;
	double xdist,ydist,zdist,xdistmax,ydistmax,zdistmax,curdist=0.;
//      Find furthest point
        for(int x=1;x<140;x++){//loop over x pixels
        	for(int y=1;y<128;y++) {//loop over y pixels
                	for(int z=1;z<timetemp;z++) {//loop over time buckets
                        	if(hcurrenttrack3D->GetBinContent(x+1,y+1,z+1)>10.) {//apply threshold
					Dist(x,y,z,xdist,ydist,zdist);
					Dist(maxx,maxy,maxz,xdistmax,ydistmax,zdistmax);
                                        curdist=(xdist-xdistmax)*(xdist-xdistmax)+(ydist-ydistmax)*(ydist-ydistmax)+(zdist-zdistmax)*(zdist-zdistmax);//distance squared from 'centre' to the point
                                        hallhits2->Fill(xdist,ydist);
                                        if(curdist>maxdist) {//new furthest point
                                		distx=x-1;
                                                disty=y-1;
                                                distz=z-1;
                                                maxdist=curdist;
                                        }
                                }
                        }
        	}
	}
        if(maxdist>0){
//		hmaxtrack->Fill(sqrt(maxdist));//plot this distance
        }
	return;
}

void Analysis::GetDot(double xposmax,double yposmax ,double zposmax ,double xdistmax ,double ydistmax,double zdistmax ) {
///////////	Do the dot product between the center of the decay and the end of the arm and every other point
		if(xposmax<=0 || yposmax<=0 || zposmax<=0 || xdistmax<=0 || ydistmax<=0 || zdistmax<=0) return;
                double v1x=0.,v2x=0.,v1y=0.,v2y=0.,v1z=0.,v2z=0.,v1t=0.,v2t=0.,dot_product=0.;
                v1x=(xposmax-xdistmax);//position (mm)
                v1y=(yposmax-ydistmax);
                v1z=(zposmax-zdistmax);
                v1t=sqrt(v1x*v1x+v1y*v1y+v1z*v1z);
                for(int x=1;x<140;x++) {
                        for(int y=1;y<128;y++) {
                                for(int z=1;z<timetemp;z++) {
                                        if(hcurrenttrack3D->GetBinContent(x+1,y+1,z+1)>0.) {
                                                double posx=0.,posy=0.,posz=0.;
						Dist(x,y,z,posx,posy,posz);
                                                v2x=posx-xdistmax;
                                                v2y=posy-ydistmax;
                                                v2z=posz-zdistmax;
                                                v2t=sqrt(v2x*v2x+v2y*v2y+v2z*v2z);
                                                dot_product=(v1x*v2x+v1y*v2y+v1z*v2z)/(v1t*v2t);
                                                double angle_sh=((4.+acos(dot_product)/d2r));
                                                if(angle_sh>180) angle_sh=360-angle_sh;//if overflowed then go back to zero
                                                if(good_evt<draw_events) hdot[good_evt]->Fill(angle_sh);
                                                hangle->Fill(angle_sh);//fill for each event then destroy (artifically shift by 4 degrees to allow zero degree peak)
						if(good_evt<draw_events) {
							double thistheta,thispsi;
							thistheta=acos(v2z/v2t);
							thispsi=atan2(v2y,v2x);
//							if(thispsi<0) thispsi+=(8.*atan(1.));//add pi
							hangle_arms[good_evt]->Fill(thistheta/d2r,thispsi/d2r);//0 -> 360
							//cout<<"point\t"<<x<<"\t"<<y<<"\t"<<z<<"\tVector"<<v2x<<"\t"<<v2y<<"\t"<<v2z<<"\tAngle\t"<<thistheta/d2r<<"\t"<<thispsi/d2r<<endl;
						}
                                        }
                                }
                        }
                }


	return;
}

double * Analysis::Get23Arm(double arm[3][2], int pxposmax, int pyposmax, int pzposmax, int pxdistmax, int pydistmax, int pzdistmax,bool doubleup) {
////////////// Find the path length for the 2nd and 3rd arm
	static double armxyz[3][3]={0.};
	static double test[9];
	double xposmax,yposmax,zposmax,xdistmax,ydistmax,zdistmax;
	Dist(pxposmax,pyposmax,pzposmax,xposmax,yposmax,zposmax);//center
	Dist(pxdistmax,pydistmax,pzdistmax,xdistmax,ydistmax,zdistmax);//arm
	int doubled_1=-1,doubled_2=-1;
	if(doubleup) {
		for(int i=0;i<3;i++) {
			for(int j=0;j<3;j++) {
				if(arm[i][0]==arm[j][0] && i<j) {
					doubled_1=i;
					doubled_2=j;
				}
			}
		}
	}
	for(int l=0;l<3;l++) {
        		if(arm[l][0]==0) {//fill long arm
				armxyz[l][0]=pxdistmax;
				armxyz[l][1]=pydistmax;
				armxyz[l][2]=pzdistmax;
                        }
                        if(arm[l][0]!=0) {//avoid long arm
                        	double fitmax=0,v1x,v1y,v1z,v1t;
                                v1x=-(xposmax-xdistmax);//position (mm)
                                v1y=-(yposmax-ydistmax);
                                v1z=-(zposmax-zdistmax);
                                v1t=sqrt(v1x*v1x+v1y*v1y+v1z*v1z);
                                for(int x=1;x<140;x++) {
                                	for(int y=1;y<128;y++) {
                                        	for(int z=1;z<timetemp;z++) {
                                                	if(hcurrenttrack3D->GetBinContent(x+1,y+1,z+1)>0.) {
                                                        	double xd,yd,zd,fit,fit2,angle,dist,dot;
                                                                Dist(x,y,z,xd,yd,zd);
//                                                              Calculate the furthest point away in this arm (ensure in arm by weighting the angle from dot-product)
                                                                dist=sqrt((xposmax-xd)*(xposmax-xd)+(yposmax-yd)*(yposmax-yd)+(zposmax-zd)*(zposmax-zd));
                                                                dot=(-(xposmax-xd)*v1x-(yposmax-yd)*v1y-(zposmax-zd)*v1z)/(dist*v1t);
                                                                angle=acos(dot)/d2r;
                                                                fit=dist*TMath::Exp(-(angle-arm[l][0])*(angle-arm[l][0])/(2.*10.*10.));///sigma=10
								fit2=1;//avoid getting the same arm again!
								if(doubleup && (l==doubled_2)) {
									double xold,yold,zold;
									Dist(armxyz[doubled_1][0],armxyz[doubled_1][1],armxyz[doubled_1][2],xold,yold,zold);
									double dist_2=(xd-xold)*(xd-xold)+(yd-yold)*(yd-yold)+(zd-zold)*(zd-zold);
									fit2=1./TMath::Exp(-(dist_2)/(2.*100.));;//avoid getting the same arm again! Negative weight to previous end
									fit=fit*fit2;
								}
                                                                if(fit>fitmax) {
                                                                        fitmax=fit;
                                                                        armxyz[l][0]=x-1;
                                                                        armxyz[l][1]=y-1;
                                                                        armxyz[l][2]=z-1;
                                                                }
                                                        }
                                                }
                                        }
                                }
                        }
		}
	for(int i=0;i<9;i++) {
		test[i]=armxyz[(i-(i%3))/3][i%3];
	}
	return test;
}


double Analysis::GetRange(double energy) {
	double range;
	range=energy;
	if(pressure==100) range/=2.;
	return range;
}
double Analysis::GetRange12C(double energy) {
	double range;
	range=2.4509*energy*energy*energy-13.578*energy*energy+37.655*energy+2.5909;
	if(pressure==100) range/=2.;
	return range;
}
double Analysis::GetRange16O(double energy) {
	double range;
	range=2.163*energy*energy*energy-12.482*energy*energy+34.396*energy+1.5331;
	if(pressure==100) range/=2.;
	return range;
}

double Analysis::GetEnergy(double range) {//alpha particles in 50 Torr CO2
	double energy;
	double a,b,c,d;
	if(pressure==100) range*=2.;
	a=-1.3733e-6;
	b=2.4464e-4;
	c=6.6056e-3;
	d=-2.8223e-2;
	energy=a*range*range*range+b*range*range+c*range+d;//MeV
	if(energy<0) energy=0.;
	return energy;
}

double Analysis::GetEnergy12N(double range) {//12N energy from range in mm
        double energy;
        double a,b,c,d;
        a=-4.212e-8;
        b=2.4587e-5;
        c=3.212e-2;
        d=-1.136;
        energy=a*range*range*range+b*range*range+c*range+d;//MeV
        return energy;
}


int Analysis::CleanY(int &maxx, int &maxy, int &maxz, int &distx, int &disty, int &distz,int beamx, int beamy, int beamz,bool &forward_arm) {
	const int n_tries=3;
	int tries_center=0,tries=0;
	bool moved_center=false;
	bool good_arm=true;
	int npeaks_ang_=0;
	RemoveOutliers(maxx,maxy,maxz);
	while(tries_center<n_tries && tries<n_tries && good_arm) {
		GetArm(maxx,maxy,maxz,distx,disty,distz);//dist(xyz) is the pixel of end of arm
		double x1,x2,y1,y2,z1,z2;
		double arm[3][2]={0.};
		bool doubleup,forward_arm;
	        Dist(maxx,maxy,maxz,x1,y1,z1);//xyz1 is position of center
	        Dist(distx,disty,distz,x2,y2,z2);///xyz2 is position of end of arm
	        double theta,psi,v1,v2,v3,vt;
	        v1=(x2-x1);//vector along arm
	        v2=(y2-y1);//vector along arm
	        v3=(z2-z1);//vector along arm
	        vt=sqrt(v1*v1+v2*v2+v3*v3);
	        theta=acos(v3/vt)/d2r;//Calculate theta psi of vector (v1,v2,v3)
	        psi=atan2(v2,v1)/d2r;
	        hthetapsi->Fill(theta,psi);
	        maxdist=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);//arm distance
	        GetDot(x2,y2,z2,x1,y1,z1);
	        for(int l=0;l<45;l++) if(hangle->GetBinContent(l+1)<=1.) hangle->SetBinContent(l+1,0);//remove low counts or peak finder will pick these out
////////////////////////////////////////////////////////////////////////////////////////////////
///		Peaks should appear in the dotproduct corresponding to the other arms which are at a fixed angle relative to the first arm vector
////////////////////////////////////////////////////////////////////////////////////////////////
	        TSpectrum *s3 = new TSpectrum();
	        const int npeaks_ang = s3->Search(hangle,2,"nobackground",0.05);//search for peaks
	        double *xpeaks_dot = s3->GetPositionX();
	        double *xpeaks_dot_h = s3->GetPositionY();
	        if(good_evt<draw_events) hdot[good_evt]->SetTitle(Form("N peaks:%i",npeaks_ang));
	        doubleup = (npeaks_ang==2) ? true : false;//if doubling up
	        for(int m=0;m<npeaks_ang;m++) {
	        	arm[m][0]=xpeaks_dot[m];
	                arm[m][1]=xpeaks_dot_h[m];
	                if(arm[m][0]<45. && arm[m][1]>5.) forward_arm=true;//If we have a forward peak with more than 5 in the peak then that is good
	        }
	        if(!forward_arm) {//bad end-point?
	        	hcurrenttrack3D->SetBinContent(distx+1,disty+1,distz+1,-hcurrenttrack3D->GetBinContent(distx+1,disty+1,distz+1));//set 'fake' max to negative (to ignore but not to lose the value) and re-do // 
	                tries++;
	                for(int m=0;m<45;m++) {
	                        hangle->SetBinContent(m,0);
	                        if(good_evt<draw_events) hdot[good_evt]->SetBinContent(m,0);
	                }
	        }
	        else {//Good end-point --> plot dist
	                good_arm=false;
	                if(maxdist>0){
	                        hmaxtrack->Fill(sqrt(maxdist));//plot this distance
	                        double height_trav=abs(z1-z2);
	                        hmaxtrack_height->Fill(sqrt(maxdist),height_trav);
	                        hlengthvsE->Fill(sqrt(maxdist),totalEnergy_glob*E_scaling+7.272);
	                }
	                npeaks_ang_=npeaks_ang;
	        }
/*
	        if(tries==n_tries) {//get rid of extraneous outliers didn't work, try changing the center
//		First add back the erroneously removed outliers
	                for(int x=0;x<140;x++) {
	                	for(int y=0;y<140;y++) {
	                        	for(int z=0;z<timetemp;z++) {
	                                		if(hcurrenttrack3D->GetBinContent(x,y,z)<0) hcurrenttrack3D->SetBinContent(x,y,z,abs(hcurrenttrack3D->GetBinContent(x,y,z)));//negative values were removed, make them positive$
	                                }
	                         }
////////////////////////////////////////////////
	                         hcurrenttrack3D->SetBinContent(maxx+1,maxy+1,maxz+1,0.25*hcurrenttrack3D->GetBinContent(maxx+1,maxy+1,maxz+1));//temporarily remove the old center by dividing it by 4
	                         GetCenter(beamx,beamy,beamz,maxx,maxy,maxz,5.0-0.5*tries_center);//each attempt, try a looser diffuseness to Gaussian, don't be so tied to the beam
	                         tries=0;
	                         tries_center++;
	                         moved_center=true;
	//                       Try once more with the moved center
	                }
	                if(distx==-1 && disty==-1 && distz==-1) return 0;
	                if(doubleup && (arm[0][0]>120. || arm[1][0]>120.)) {//overlapping??
	                        arm[2][0]=max(arm[0][0],arm[1][0]);//take highest angle
	                        arm[2][1]=min(arm[0][1],arm[1][1]);//take smaller peak (not long arm)
	                }
	        }
*/
	}
	return npeaks_ang_;
}


void Analysis::RemoveOutliers(int maxx,int maxy,int maxz) {
	TH1F *hdistances = new TH1F("hdistances","Distance to points",100,0,200.);
	double distx_b,disty_b,distz_b,distx_p,disty_p,distz_p;
	for(int x=0;x<140;x++) {
		for(int y=0;y<128;y++) {
			for(int z=0;z<timetemp;z++) {
				if(hcurrenttrack3D->GetBinContent(x+1,y+1,z+1)>10.) {
					double distance_to_point=0;
					Dist(maxx,maxy,maxz,distx_b,disty_b,distz_b);
					Dist(x,y,z,distx_p,disty_p,distz_p);
					distance_to_point=sqrt((distx_b-distx_p)*(distx_b-distx_p)+(disty_b-disty_p)*(disty_b-disty_p)+(distz_b-distz_p)*(distz_b-distz_p));
					hdistances->Fill(distance_to_point);//get histogram of points
				}
			}
		}
	}
	//remove points based on hdistances distribution
//	for(int i=0;i<100;i++) cout<<i<<"\t"<<hdistances->GetBinContent(i+1)<<endl;
	double std_dev = hdistances->GetStdDev(1);//std dev for x axis
	double mean = hdistances->GetMean(1);//mean for x axis
        for(int x=0;x<140;x++) {
                for(int y=0;y<128;y++) {
                        for(int z=0;z<timetemp;z++) {
                                if(hcurrenttrack3D->GetBinContent(x+1,y+1,z+1)>10.) {
                                        double distance_to_point=0;
                                        Dist(maxx,maxy,maxz,distx_b,disty_b,distz_b);
                                        Dist(x,y,z,distx_p,disty_p,distz_p);
                                        distance_to_point=sqrt((distx_b-distx_p)*(distx_b-distx_p)+(disty_b-disty_p)*(disty_b-disty_p)+(distz_b-distz_p)*(distz_b-distz_p));
//					cout<<x<<"\t"<<y<<"\t"<<z<<"\t"<<distance_to_point<<"\t"<<mean+2.*std_dev<<endl;
					if(distance_to_point>mean+2.*std_dev) {
						hcurrenttrack3D->SetBinContent(x+1,y+1,z+1,0.);//if more than 2 std. dev. away then remove
//						cout<<"REMOVE\t"<<x<<"\t"<<y<<"\t"<<z<<endl;
					}
                                }
                        }
                }
        }

	delete hdistances;
	return;
}


void Analysis::GetPlane() {//takes two vectors and gives orthogonal vector which defines the plane
//For vectors {v1x,v1y,v1z} and {v2x,v2y,v2z}, the co-orthogonal vector is given by v1 x v2
	double v1x=0.,v2x=0.,v3x=0.,v1y=0.,v2y=0.,v3y=0.,v1z=0.,v2z=0.,v3z=0.,v3t=0.;
	v1x = HoughLine[0][0];
	v1y = HoughLine[0][1];
	v1z = HoughLine[0][2];
	v2x = HoughLine[1][0];
	v2y = HoughLine[1][1];
	v2z = HoughLine[1][2];
	v3x = v1y*v2z-v1z*v2y;
	v3y = v1z*v2x-v1x*v2z;
	v3z = v1x*v2y-v1y*v2x;
	v3t = sqrt(v3x*v3x+v3y*v3y+v3z*v3z);
	Normal[0]=v3x/v3t;
	Normal[1]=v3y/v3t;
	Normal[2]=v3z/v3t;
//
        v1x = HoughLine[2][0];
        v1y = HoughLine[2][1];
        v1z = HoughLine[2][2];
        v2x = HoughLine[1][0];
        v2y = HoughLine[1][1];
        v2z = HoughLine[1][2];
        v3x = v1y*v2z-v1z*v2y;
        v3y = v1z*v2x-v1x*v2z;
        v3z = v1x*v2y-v1y*v2x;
        v3t = sqrt(v3x*v3x+v3y*v3y+v3z*v3z);
//
//	Normal_1.Normal_2
	v3x/=v3t;
	v3y/=v3t;
	v3z/=v3t;
	double dotproduct;
	dotproduct=Normal[0]*v3x+Normal[1]*v3y+Normal[2]*v3z;
	if(good_evt<draw_events) hdistxy[good_evt]->SetTitle(Form("Dot product %g",dotproduct));
	cout<<"theta "<<acos(v3x)/d2r<<" psi "<<atan(v3y/v3x)/d2r<<endl;
	return;
}


void Analysis::Fits() {

	fit = new TF1("fit","[0]*exp(-x/[1])",0,30);
	fit->SetParameter(0,hD2Ptime->GetBinContent(1));
	fit->SetParameter(1,11.);//t1/2=11 ms
	cout<<"Fitting half-life"<<endl;
	hD2Ptime->Fit("[0]*exp(-x/[1])");
//	hD2Ptime->Fit("fit","0");
	hD2Ptime->Draw();
	fit->Draw("SAME");
	cout<<"Fitted"<<endl;
	return;
}


void Analysis::Flatten(int maxx, int maxy, int maxz) {
//Plane defined by normal vector (A,B,C) and point Q = (x_0,y_0,z_0)
//Equation for plane is A(x-x_0) + B(y-y_0) + C(z-z_0) = 0
//Or Ax + By + Cz + D = 0 with D = -Ax_0 - By_0 - Cz_0
///////////////
//Distance from point (P) to plane is |A(x_1-x_0)+B(y_1-y_0)+C(z_1-z_0)|/sqrt(A^2+B^2+C^2)
//Or |Ax_1 + By_1 + Cz_1 + D|/sqrt(A^2+B^2+C^2)
///////////////
	int x_low=140,x_high=0;
	double Q_x1=0.,Q_y1=0.,Q_z1=0.,Q_x2=0.,Q_y2=0.,Q_z2=0.;
	double pix_x1=0.,pix_x2=0.,pix_y1=0.,pix_y2=0.,pix_z1=0.,pix_z2=0.;
	double d_x=0.,d_y=0.,d_z=0.;
//	Pin to plane via center maxx,maxy,maxz -> Q=(maxx,maxy,maxz)
//	One vector v1 = (x_high-x_low,y_high-y_low,z_high-z_low)
//	Normal vector: n = v1 X v2
	double n[3];
	n[0]=Normal[0];
	n[1]=Normal[1];
	n[2]=Normal[2];
//	Nearest point on plane (G) to point P (x,y,z)
//	Where P - lambda*n intersects plane
//	A,B,C are given by n
	double nt=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	n[0]=n[0]/nt;
	n[1]=n[1]/nt;
	n[2]=n[2]/nt;
	cout<<"Normal to plane "<<n[0]<<"\t"<<n[1]<<"\t"<<n[2]<<endl;
	double A,B,C,D,Ap,Bp,Cp,Dp;
	A=n[0];
	B=n[1];
	C=n[2];
	Ap=A;
	Bp=B;
	Cp=C;
	double pos_maxx = decay_x,pos_maxy=decay_y,pos_maxz=decay_z;
	D=-A*pos_maxx - B*pos_maxy - C*pos_maxz;
//	A(P_x-lambda*n_x) + B(P_y-lambda*n_y) + C(P_z-lambda*n_z) + D = 0
	double sumlambda=-A*A-B*B-C*C;
	double consts;
	double theta_y = 2.*atan(1.) + atan2(Cp,Ap);
	cout<<"Rotate by theta_y "<<theta_y/d2r<<endl;

//	cout<<"ABC pre: "<<A<<"\t"<<B<<"\t"<<C<<endl;
	n[0]=Ap*cos(theta_y)+Cp*sin(theta_y);
        n[1]=Bp;
        n[2]=-Ap*sin(theta_y)+Cp*cos(theta_y);
        A=n[0];
        B=n[1];
	C=n[2];
	double theta_x = atan2(B,C);
	cout<<"Rotate by theta_x "<<theta_x/d2r<<endl;
//	cout<<"ABC1: "<<A<<"\t"<<B<<"\t"<<C<<endl;

        n[0]=A;
        n[1]=B*cos(theta_x)-C*sin(theta_x);
        n[2]=B*sin(theta_x)+C*cos(theta_x);
        A=n[0];
        B=n[1];
        C=n[2];
//        cout<<"ABC2 "<<A<<"\t"<<B<<"\t"<<C<<endl;
//	cout<<"LOOP!"<<endl<<endl<<endl<<endl<<endl;
	double dist_cent_x,dist_cent_y,dist_cent_z;
	Dist(decay_x,decay_y,decay_z,dist_cent_x,dist_cent_y,dist_cent_z);
	for(int x=0;x<140;x++) {
		for(int y=0;y<128;y++) {
			for(int z=0;z<timetemp;z++) {
				if(hcurrenttrack3D->GetBinContent(x+1,y+1,z+1)>0) {
					for(int i=0;i<100;i++) {
						Dist(x,y,z,d_x,d_y,d_z);
						d_x -= dist_cent_x;
						d_y -= dist_cent_y;
						d_z -= dist_cent_z;
						consts=Ap*d_x+Bp*d_y+Cp*d_z + D;
						double lambda=consts/sumlambda;
//	New point is now P-lambda*n
//						theta_x=45*d2r;
//						theta_y=45*d2r;
						lambda=0.;
						double xnew=d_x-lambda*Ap;
						double ynew=d_y-lambda*Bp;
						double znew=d_z-lambda*Cp;
//						cout<<"X: "<<d_x<<"->"<<xnew<<" Y: "<<d_y<<"->"<<ynew<<" Z: "<<d_z<<"->"<<znew<<endl;
//	Project onto 2D plane
//	Rotate such that n_x = 0 around y-axis
//						cout<<"Rotate about y by "<<theta_y*57.1<<endl;
						double xnew2=xnew*cos(theta_y) + znew*sin(theta_y);
						double ynew2=ynew;
						double znew2=-xnew*sin(theta_y)+znew*cos(theta_y);
						n[0]=Ap*cos(theta_y)+Cp*sin(theta_y);
						n[1]=Bp;
						n[2]=-Ap*sin(theta_y)+Cp*cos(theta_y);
						A=n[0];
						B=n[1];
						C=n[2];
//						cout<<"NEW1: "<<xnew2<<"\t"<<ynew2<<"\t"<<znew<<endl;

//	Rotate such that n_y = 0 around x-axis
//						cout<<"ABC "<<A<<"\t"<<B<<"\t"<<C<<endl;
//						cout<<"Rotate by theta_y"<<theta_y/d2r<<"\t"<<xnew2<<"\t"<<ynew2<<"\t"<<znew2<<endl;
//						double theta_x = atan2(B,C);
						xnew=xnew2;
						ynew=ynew2*cos(theta_x)-znew2*sin(theta_x);
						znew=ynew2*sin(theta_x)+znew2*cos(theta_x);
//	Calculate projection angle to maintain distances
//						cout<<"Rotate by theta_x"<<theta_x/d2r<<"\t"<<xnew<<"\t"<<ynew<<"\t"<<znew<<endl;
                                                n[0]=A;
                                                n[1]=B*cos(theta_x)-C*sin(theta_x);
                                                n[2]=B*sin(theta_x)+C*cos(theta_x);
                                                A=n[0];
                                                B=n[1];
                                                C=n[2];
//						cout<<"ABC2 "<<A<<"\t"<<B<<"\t"<<C<<endl;
//						cout<<"NEW2: "<<xnew<<"\t"<<ynew<<"\t"<<znew<<endl;
//						projectionangles=1.;
//						cout<<xnew<<"\t"<<ynew<<"\t"<<znew<<endl;
						hflat[good_evt]->Fill((xnew),(ynew),hcurrenttrack3D->GetBinContent(x+1,y+1,z+1)*0.01);//Project
						hflat_all->Fill((xnew),(ynew),hcurrenttrack3D->GetBinContent(x+1,y+1,z+1)*0.01);//Project
//						cout<<xnew/projectionangles<<"\t"<<ynew/projectionangles<<"\t"<<hcurrenttrack3D->GetBinContent(x+1,y+1,z+1)<<endl;
//						cout<<projectionangles<<"\t\t"<<A<<"\t"<<B<<"\t"<<C<<endl;
					}
				}
			}
		}
	}
	return;
}


void Analysis::FindReactionPlane2D() {//Find reaction plane (psi only as beam path defines part of path)

//	Not used -- take the HoughTransform plane

	return;
}


void Analysis::NearestNeighbour() {
	for(int x=0;x<140;x++) {
		for(int y=0;y<128;y++) {
			for(int z=0;z<timetemp;z++) {
				//For each voxel find the nearest neighbour
				if(hcurrenttrack3D->GetBinContent(x+1,y+1,z+1)>0.) {
					for(int dstep=1;dstep<=20;dstep++) {//how many voxels away?
						if(dstep==20) hcurrenttrack3D->SetBinContent(x+1,y+1,z+1,0);//Remove this outlier as it is more than 20 voxels away
						for(int xtest=-dstep;xtest<=dstep;xtest++) {
							for(int ytest=-dstep;ytest<=dstep;ytest++) {
									for(int ztest=-dstep;ztest<=dstep;ztest++) {
//									cout<<"Testing distance "<<dstep<<"   "<<x<<","<<y<<","<<z<<" versus "<<x+xtest<<","<<y+ytest<<","<<z+ztest<<endl;
									bool self=false;
									if(xtest==0 && ytest==0 && ztest==0) self=true;
									if(hcurrenttrack3D->GetBinContent(x+xtest+1,y+ytest+1,z+ztest+1)>0 && !self) {
										//Found nearest neighbour
										if(good_evt<draw_events) hnearestneighbour[good_evt]->Fill(dstep);
										//Exit loops
										xtest=100;
										ytest=100;
										ztest=100;
										dstep=100;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return;
}


void Analysis::FancyCenter() {
//	Use the Hough transform parameters to find the center
//	This is the point where the main 2 lines overlap the best

	return;
}


void Analysis::GetOverlapPoint(double& x0,double& y0,double& z0) {
//	Get the point of overlap between the side region and the central region
	int selside = (canstripmatchL ? 5 : 0);
	double weighted_position=0.,total_weight=0,weighted_time=0.;
	for(int y=0;y<128;y++) {
		if(centralhitsevt[selside][y]>0.) {//adjoining hit
			weighted_position+=y*centralhitsevt[selside][y];
			weighted_time+=centralhitsTevt[selside][y]*centralhitsevt[selside][y];
			total_weight+=centralhitsevt[selside][y];
		}
	}
	cout<<weighted_position<<"\t"<<total_weight<<endl;
	weighted_position/=total_weight;
	cout<<"Weighted position "<<weighted_position<<endl;
	x0=selside+64;
	y0=weighted_position;
	z0=weighted_time/total_weight;
	return;
}


void Analysis::StripsChains() {
		int bestx=-1,besty=-1,bestl=-1,bestm=-1;
		double besttimediff=50.,energy=0.,bestt=0.;
////////	Test to see how flat the strips and chains are in time, if they are too flat then try use different method
///		Separate into beam left and beam right
		hchainscheckL = new TH1F("hchainscheckL","Time of chains",timetemp,0,timetemp);
		hchainscheckR = new TH1F("hchainscheckR","Time of chains",timetemp,0,timetemp);
		for(int x=0;x<striparsize;x++) {
			for(int m=0;m<2;m++) {
				if(mmchainsT[x][m]>105.) {
					if(x<=66. )hchainscheckL->Fill(mmchainsT[x][m]);
					if(x>66. ) hchainscheckR->Fill(mmchainsT[x][m]);
					if(good_evt<draw_events) hchainsevt[good_evt]->Fill(mmchainsT[x][m]);
				}
			}
		}
		canstripmatchL=true;
		canstripmatchR=true;
		if(hchainscheckL->GetMaximum()<8. || fulleventreadout) {//can do the normal method
			canstripmatchL=true;
		}
		if(hchainscheckR->GetMaximum()<8. || fulleventreadout) {//can do the normal method
			canstripmatchR=true;
		}
		delete hchainscheckL;
		delete hchainscheckR;
		if(canstripmatchL||canstripmatchR) {//check EITHER side is working
////////////FILL best hit for each chain
 			int side=0;
			maxdist=0.;
			for(int y=0;y<striparsize;y++) {//loop over strip
				for(int m=0;m<2;m++){//loop over hit for strip
					besttimediff=10.;
					bestx=0;
					besty=0;
					for(int x=0;x<striparsize;x++) {
						side = (x<66) ? 0 : 1;//if x is less than 66  then left side
						if((side==0 && canstripmatchL) || (side==1 && canstripmatchR)) {//make sure we are matching strips and chains on a side that will work
							for(int l=0;l<2;l++) {
								if(mmchainsE[x][l]>0. && mmstripsE[y][side][m]>0. && mmchainsT[x][l]>0 && mmstripsT[y][side][m]>0) {//only bother for non-zero elements
									if(abs(mmchainsT[x][l]-mmstripsT[y][side][m])<besttimediff) {//new closer matching
											besttimediff=abs(mmchainsT[x][l]-mmstripsT[y][side][m]);
											bestx=x;
											besty=y;
											bestl=l;
											bestm=m;
//											energy=0.5*(mmstripsE[x][l]+mmchainsE[y][m]);
											energy=max(mmchainsE[x][l],mmstripsE[y][side][m]);
											bestt=(mmstripsT[y][side][m]);
//											bestt=0.5*(mmchainsT[x][l]+mmstripsT[y][m]);
											if(besttimediff<1.){//If it's a good match then try drawing anyway
//												if(good_evt<draw_events) htrackxy[good_evt]->Fill(bestx,besty,energy);
//												if(good_evt<draw_events) htrackyz[good_evt]->Fill(besty,bestt,energy);
//												if(good_evt<draw_events) htrackxz[good_evt]->Fill(bestx,bestt,energy);
//												if(good_evt<draw_events) htrack3D[good_evt]->Fill(bestx,besty,bestt,energy);
//												if(good_evt<draw_events) hsuper[good_evt]->Fill(bestx,besty);
//												hcurrenttrack3D->Fill(bestx,besty,bestt,energy);
//												hcurrenttrack3D_2->Fill(bestx,besty,bestt,min(mmchainsE_2[x][l],mmstripsE_2[y][side][m]));
//												hcurrentYZ->Fill(besty,bestt,energy);
											}
									}
								}
							}
						}
					}
					if(bestx>0 && besty>0) {//found a suitable match
						if(good_evt<draw_events) htrackxy[good_evt]->Fill(bestx,besty,energy);
						if(good_evt<draw_events) hsuper[good_evt]->Fill(bestx,besty);
						if(good_evt<draw_events) hdeltat[good_evt]->Fill(bestx,besty,besttimediff);
						if(good_evt<draw_events) htrackyz[good_evt]->Fill(besty,bestt,energy);
						if(good_evt<draw_events) htrackxz[good_evt]->Fill(bestx,bestt,energy);
//						if(good_evt<draw_events) htrack3D[good_evt]->Fill(bestx,besty,bestt,energy);
						energy=max(mmchainsE[bestx][bestl],mmstripsE[besty][side][bestm]);
						hcurrenttrack3D->Fill(bestx,besty,bestt,energy);

						hcurrentYZ->Fill(besty,bestt,energy);

//						hcurrenttrack3D_2->Fill(bestx,besty,bestt,min(mmchainsE_2[bestx][bestl],mmstripsE_2[besty][side][bestm]));
//						if(good_evt<draw_events) htrackxysmooth[good_evt]->Fill(bestx,besty,min(mmchainsE_2[bestx][bestl],mmstripsE_2[besty][side][bestm]));
						besttimediff=10.;
						bestx=-1;
						besty=-1;
					}
				}
			}
////////////FILL best hit for each chain
                	for(int x=0;x<striparsize;x++) {
                		for(int l=0;l<2;l++) {
                	                besttimediff=10.;
					bestx=0;
					besty=0;
					side = (x<66) ? 0 : 1;//if x is less than 66  then left side
					if((side==0 && canstripmatchL) || (side==1 && canstripmatchR)) {//make sure we are matching strips and chains on a side that will work
                				for(int y=0;y<striparsize;y++) {
                	        			for(int m=0;m<2;m++){
                	                	                if(mmchainsE[x][l]>0. && mmstripsE[y][side][m]>0. && mmchainsT[x][l]>0 && mmstripsT[y][side][m]>0) {//only bother for non-zero elements
                	                	                         if(abs(mmchainsT[x][l]-mmstripsT[y][side][m])<besttimediff) {//new closer matching
                	                	                                        besttimediff=abs(mmchainsT[x][l]-mmstripsT[y][side][m]);
                	                	                                        bestx=x;
                	                	                                        besty=y;
											bestl=l;
											bestm=m;
//                     	                	                                   energy=0.5*(mmstripsE[x][l]+mmchainsE[y][m]);
                        	        	                                        energy=max(mmchainsE[x][l],mmstripsE[y][side][m]);
											bestt=(mmchainsT[x][l]);
	                                                                        }

                        	        	                                        if(fabs(mmchainsT[x][l]-mmstripsT[y][side][m])<=1.){//If it's a good match then try drawing anyway
  //                      	        	                                                if(good_evt<draw_events) htrackxy[good_evt]->Fill(bestx,besty,energy);
    //                    	        	                                                if(good_evt<draw_events) hsuper[good_evt]->Fill(bestx,besty);
      //                  	        	                                                if(good_evt<draw_events) htrackyz[good_evt]->Fill(besty,bestt,energy);
        //                	        	                                                if(good_evt<draw_events) htrackxz[good_evt]->Fill(bestx,bestt,energy);
	//	                      	               	                                        if(good_evt<draw_events) htrack3D[good_evt]->Fill(bestx,besty,bestt,energy);
          //              	                	                                        hcurrenttrack3D->Fill(bestx,besty,bestt,energy);
	//											hcurrentYZ->Fill(besty,bestt,energy);

											}
	 	                                                       }
	 	                                               }
	 	                                       }
	                                }
	                                if(bestx>0 && besty>0) {//found a suitable match
	                                        if(good_evt<draw_events) htrackxy[good_evt]->Fill(bestx,besty,energy);
	                                        if(good_evt<draw_events) hsuper[good_evt]->Fill(bestx,besty);
						if(good_evt<draw_events) htrackyz[good_evt]->Fill(besty,bestt,energy);
						if(good_evt<draw_events) htrackxz[good_evt]->Fill(bestx,bestt,energy);
//						if(good_evt<draw_events) htrack3D[good_evt]->Fill(bestx,besty,bestt,energy);
						energy=max(mmchainsE[bestx][bestl],mmstripsE[besty][side][bestm]);
						hcurrenttrack3D->Fill(bestx,besty,bestt,energy);
						hcurrentYZ->Fill(besty,bestt,energy);

						hcurrenttrack3D_2->Fill(bestx,besty,bestt,min(mmchainsE_2[bestx][bestl],mmstripsE_2[besty][side][bestm]));
						if(good_evt<draw_events) htrackxysmooth[good_evt]->Fill(bestx,besty,min(mmchainsE_2[bestx][bestl],mmstripsE_2[besty][side][bestm]));
	                                        bestx=0;
	                                        besty=0;
	                                }
	                        }
	                }
		}
		if(!canstripmatchL || !canstripmatchR) {//Other method for strip/chain matching needed as time differences are too small!
			cout<<evtno<<"\tOther method activate!"<<canstripmatchL<<"\t"<<canstripmatchR<<endl;
			if(canstripmatchL && canstripmatchR) cout<<BOLDRED<<"Both sides are bad -- check this event"<<RESET<<endl;
			int smallestch=200,largestch=0,smallestst=200,largestst=0;//Largest/smallest chain/strip
			for(int l=0;l<striparsize;l++) {//find range of strips
				for(int m=0;m<2;m++) {
					int side = (canstripmatchL ? 1 : 0);//stripmatch L --> side = 0
					if(mmstripsT[l][side][m]>105.) {//Good strip event
						if(l<smallestst) smallestst=l;
						if(l>largestst) largestst=l;
					}
				}
			}
			for(int l=0;l<striparsize;l++) {//find range of chains
				for(int m=0;m<2;m++) {
					bool correctside=false;
					if((canstripmatchL && l>66) || (canstripmatchR && l<66)) correctside=true;
					if(mmchainsT[l][m]>105. && correctside) {//Good strip event
						if(l<smallestch) smallestch=l;
						if(l>largestch) largestch=l;
					}
				}
			}
			cout<<evtno<<"\t"<<smallestch<<"\t"<<largestch<<"\t"<<smallestst<<"\t"<<largestst<<endl;
			bool choice = ((largestch-smallestch) > (largestst-smallestst))? true : false;//if the range of chains is higher than the range of strips then step over these to give a smoother line
			double x0,y0,x1,y1,z0;//defines a point on the line - given by the continuous overlap with the central pad
			GetOverlapPoint(x0,y0,z0);
			if(fabs(x0-smallestch)<(x0-largestch)) {
				x1=largestch;//if largestx is further away from the center then this is the endpoint
			}
			else {
				x1=smallestch;
			}
			if(fabs(y0-smallestst)<(y0-largestst)) {
				y1=largestst;
			}
			else {
				y1=smallestst;
			}
			cout<<"(x0,y0) = "<<x0<<","<<y0<<"\t(x1,y1) = "<<x1<<","<<y1<<" time = "<<z0<<endl;
			double slope = (y1-y0)/(x1-x0);
			double c= y0-slope*x0;
			cout<<"slope: "<<slope<<" c: "<<c<<endl;

			if(choice) {//loop over the chains
				cout<<"Draw chains"<<endl;
				for(int l=0;l<striparsize;l++) {
					for(int m=0;m<2;m++) {
	                                        bool correctside=false;
               		                        if((canstripmatchL && l>66) || (canstripmatchR && l<66)) correctside=true;

						if(mmchainsT[l][m]>50. && correctside) {
							double chains_x = 0;
//							y=mx+c -> x=(y-c)/,
//							chains_x = (l-c)/slope;
							chains_x = slope*l+c;
							if(good_evt<draw_events) htrackxy[good_evt]->Fill(l,chains_x,mmchainsE[l][m]);
							hcurrenttrack3D->Fill(l,chains_x,mmchainsT[l][m],mmchainsE[l][m]);
						}
					}
				}
			}
			else{//loop over the strips instead
				cout<<"Draw strips"<<endl;
                                for(int l=0;l<striparsize;l++) {
                                        for(int m=0;m<2;m++) {
                                                int side=0;
						if(l>66) side=1;
                                                if(mmstripsT[l][side][m]>50.) {
                                                        double strips_y = 0;
//                                                      y=mx+c -> x=(y-c)/,
//                                                        strips_y = slope*l + c;
                                                        strips_y = (l-c)/slope;
                                                        if(good_evt<draw_events) htrackxy[good_evt]->Fill(strips_y,l,mmstripsE[l][side][m]);
							hcurrenttrack3D->Fill(strips_y,l,mmstripsT[l][side][m],mmstripsE[l][side][m]);

                                                }
                                        }
                                }

			}
		}
/////////////
		return;
}



void Analysis::OtherArms() {
/*
		double armxyz[3][3]={0.};
		double ArmLength[3]={0.};
		if(debug) cout<<"Do some peaky stuff"<<endl;
		if((npeaks_ang_==3 || npeaks_ang_==2) && forward_arm) {//everything is OK if 3 peaks and the first is below 45
//			Check the peaks are big enough (and proper peaks)
			bool big_enough=true;
			for(int l=0;l<3;l++) {
				if(arm[l][1]<3.) {
					big_enough=false;
				}
			}
			if(big_enough) {
//				Shift to match first
				double shift_val=min(arm[0][0],min(arm[1][0],arm[2][0]));
				for(int l=0;l<3;l++) {
					arm[l][0]-=shift_val;
				}
				for(int l=0;l<3;l++) {
					for(int m=0;m<3;m++) {
						for(int k=0;k<3;k++) {
							if(k!=l && m!=l && arm[k][0]<arm[l][0] && arm[k][0]<arm[m][0]) 	{//don't self match and k is the smallest
								htheta2theta3->Fill(arm[l][0],arm[m][0]);
								htheta2theta3->Fill(arm[m][0],arm[l][0]);
							}
						}
					}
				}
				if(debug) cout<<"Get 23 arm"<<endl;
				double *armxyz_ = Get23Arm(arm,maxx,maxy,maxz,distx,disty,distz,doubleup);//return 1D array
				if(debug) cout<<"Got 2nd 3rd arms"<<endl;
				goodarms=true;
				for(int l=0;l<3;l++) {//copy to 2D array
					for(int m=0;m<3;m++) {
						armxyz[l][m]=*(armxyz_+3*l+m);
					}
				}
				for(int l=0;l<3;l++) {
					double dist_=0,x_1,x_2,y_1,y_2,z_1,z_2;
					Dist(armxyz[l][0],armxyz[l][1],armxyz[l][2],x_1,y_1,z_1);//arm pixel -> arm position
					Dist(maxx,maxy,maxz,x_2,y_2,z_2);//center pixel -> arm position
					dist_=sqrt((x_1-x_2)*(x_1-x_2)+(y_1-y_2)*(y_1-y_2)+(z_1-z_2)*(z_1-z_2));
					if(debug)cout<<i<<"\t"<<armxyz[l][0]<<"\t"<<armxyz[l][1]<<"\t"<<armxyz[l][2]<<"\t"<<maxx<<"\t"<<maxy<<"\t"<<maxz<<"\t"<<dist_<<endl;
					hlength_all->Fill(dist_);
					ArmLength[l]=dist_;
				}
			}
		}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(goodarms || 1) {
			double totallength=0,E_a[3],eps_a[3],Etotal=0.;
			for(int l=0;l<3;l++) {
				totallength+=ArmLength[l];
				E_a[l]=GetEnergy(ArmLength[l]);
				Etotal+=E_a[l];
			}
			hEx->Fill(Etotal+7.272);
			for(int l=0;l<3;l++) eps_a[l]=E_a[l]/Etotal;
			double x,y;
//			TODO: randomly order
			x=(eps_a[1]-eps_a[0])/sqrt(3.);
			y=(2.*eps_a[2]-eps_a[1]-eps_a[0])/3.;
			hDalitz->Fill(x,y);
			htotal_length->Fill(totallength);
		}
*/
	return;
}


void Analysis::DefineEventHistograms() {
		ostringstream name_h,name2_h,name_side,name_E,name_pulse_chn,name_pulse_chndc,name_track,name_trackxy,name_trackyz,name_beam,name_nn,name_mult,name_projXY,name_projYZ,name_projXZ;
		ostringstream name_dot,name_trackxz,name_dedx,name_angarms,name_super,name_waveform,name_delta,name_energy,name_flat,name_arm1_xy,name_arm1_yz;//variables for naming hists
		ostringstream name_HoughXY,name_HoughYZ,name_HoughXZ,name_fit_beamXY,name_fit_beamYZ,name_HoughXY2,name_HoughYZ2,name_HoughXZ2,name_distxy,name_distxz,name_distyz,name_stage1,name_stage2,name_stage3,name_chains,name_waveform_hit,name_maxchains,name_maxstrips;
		if(good_evt<draw_events) {//define histograms which are done for each good event (up until draw_events)
			if(hpulse[good_evt]!=0 || henergy[good_evt]!=0) {//delete old histograms first
				delete hpulse[good_evt];
				delete hwaveform[good_evt];
				delete hwaveform_dc[good_evt];
				delete hwaveform_hit[good_evt];
				delete hpulse_b[good_evt];
				delete hHits[good_evt];
				delete hHits_side[good_evt];
				delete hpulse_chn[good_evt];
				delete hpulse_chndc[good_evt];
//				delete htrack3D[good_evt];
				delete henergy[good_evt];
				delete htrackxy[good_evt];
				delete htrackyz[good_evt];
				delete htrackxz[good_evt];
				delete hbeam[good_evt];
				delete hbeam2[good_evt];
				delete hfitbeamXY[good_evt];
				delete hfitbeamYZ[good_evt];
				delete htrackxysmooth[good_evt];
				delete hdot[good_evt];
				delete hdedx[good_evt];
				delete hangle_arms[good_evt];
				delete hsuper[good_evt];
				delete hdeltat[good_evt];
				delete hflat[good_evt];
				delete hnearestneighbour[good_evt];
				delete harm1_xy[good_evt];
				delete harm1_yz[good_evt];
				delete hmults[good_evt];
				delete hHoughXY[good_evt];
				delete hHoughXZ[good_evt];
				delete hHoughYZ[good_evt];
				delete hHoughXY2[good_evt];
				delete hHoughXZ2[good_evt];
				delete hdistxy[good_evt];
				delete hdistyz[good_evt];
				delete hdistxz[good_evt];
				delete hdistxy1[good_evt];
				delete hdistyz1[good_evt];
				delete hdistxz1[good_evt];
				delete hdistxy2[good_evt];
				delete hdistyz2[good_evt];
				delete hdistxz2[good_evt];
				delete hdistxy2a[good_evt];
				delete hdistyz2a[good_evt];
				delete hdistxz2a[good_evt];
				delete hdistxy2b[good_evt];
				delete hdistyz2b[good_evt];
				delete hdistxz2b[good_evt];
				delete hdistxy2c[good_evt];
				delete hdistyz2c[good_evt];
				delete hdistxz2c[good_evt];
				delete hdistxy3[good_evt];
				delete hdistyz3[good_evt];
				delete hdistxz3[good_evt];
				delete harmxy[good_evt];
				delete harmyz[good_evt];
				delete harmxz[good_evt];
				delete hHoughprojXY[good_evt];
				delete hHoughprojYZ[good_evt];
				delete hHoughprojXZ[good_evt];
				delete hstage1[good_evt];
				delete hstage2[good_evt];
				delete hstage3[good_evt];
				delete hchainsevt[good_evt];
				delete hmaxstrips[good_evt];
				delete hmaxchains[good_evt];
			}
			name2_h<<"hpulse"<<evtno;
			hpulse[good_evt]=new TH2F(name2_h.str().c_str(),name2_h.str().c_str(),timetemp,0,timetemp,420,-200,4000);
			hpulse[good_evt]->SetDirectory(0);
			name2_h<<"_b";
			hpulse_b[good_evt]=new TH2F(name2_h.str().c_str(),name2_h.str().c_str(),timetemp,0,timetemp,400,0,4000);
			hpulse_b[good_evt]->SetDirectory(0);
			name_h<<"hHits"<<evtno;
			hHits[good_evt]=new TH2F(name_h.str().c_str(),name_h.str().c_str(),140,0,140,128,0,128);
			name_side<<"hHits_side"<<evtno;
			hHits_side[good_evt]=new TH2F(name_side.str().c_str(),name_side.str().c_str(),200,0,200,200,0,200);
			name_pulse_chn<<"hpulse_chn"<<evtno;
			hpulse_chn[good_evt]=new TH1F(name_pulse_chn.str().c_str(),name_pulse_chn.str().c_str(),timetemp,0,timetemp);
			name_pulse_chndc<<"hpulse_chn_dc"<<evtno;
			hpulse_chndc[good_evt]=new TH1F(name_pulse_chndc.str().c_str(),name_pulse_chndc.str().c_str(),timetemp,0,timetemp);
			hpulse_chndc[good_evt]->SetDirectory(0);
			name_track<<"htrack3D"<<evtno;
//			htrack3D[good_evt] = new TH3F(name_track.str().c_str(),name_track.str().c_str(),140*0.5,0,140,128*0.5,0,128,128,0,timetemp);
//			htrack3D[good_evt]->SetDirectory(0);
			name_trackxy<<"htrackxy"<<evtno;
			htrackxy[good_evt] = new TH2F(name_trackxy.str().c_str(),name_trackxy.str().c_str(),140,0,140,128,0,128);
			htrackxy[good_evt]->SetDirectory(0);
			name_trackxy<<"smooth";
			htrackxysmooth[good_evt] = new TH2F(name_trackxy.str().c_str(),name_trackxy.str().c_str(),140,0,140.,128.,0,128.);
			htrackxysmooth[good_evt]->SetDirectory(0);
			name_trackyz<<"htrackyz"<<evtno;
			htrackyz[good_evt] = new TH2F(name_trackyz.str().c_str(),name_trackyz.str().c_str(),140,0,140,timetemp,0,timetemp);
			htrackyz[good_evt]->SetDirectory(0);
			name_trackxz<<"htrackxz"<<evtno;
			htrackxz[good_evt] = new TH2F(name_trackxz.str().c_str(),name_trackxz.str().c_str(),140,0,140,timetemp,0,timetemp);
			htrackxz[good_evt]->SetDirectory(0);
			name_beam<<"hbeam"<<evtno;
			hbeam[good_evt] = new TH2F(name_beam.str().c_str(),name_beam.str().c_str(),140,0,140,128,0,128);
			hbeam[good_evt]->SetDirectory(0);
			name_beam<<"_2";
			hbeam2[good_evt] = new TH2F(name_beam.str().c_str(),name_beam.str().c_str(),128,0,128,256,0,timetemp);
			name_fit_beamXY<<"hfitbeamXY"<<evtno;
			hfitbeamXY[good_evt] = new TH2F(name_fit_beamXY.str().c_str(),name_fit_beamXY.str().c_str(),140,-245,245,256,-224,224);
			name_fit_beamYZ<<"hfitbeamYZ"<<evtno;
			hfitbeamYZ[good_evt] = new TH2F(name_fit_beamYZ.str().c_str(),name_fit_beamYZ.str().c_str(),140,-140,140,timetemp,0,timetemp);
			name_dot<<"hdot"<<evtno;
			hdot[good_evt] = new TH1F(name_dot.str().c_str(),name_dot.str().c_str(),45,0,180.);
			hdot[good_evt]->SetDirectory(0);
			name_dedx<<"hdedx"<<evtno;
			hdedx[good_evt] = new TH1F(name_dedx.str().c_str(),name_dedx.str().c_str(),npoints_dedx,0,100);
			hdedx[good_evt]->SetDirectory(0);
			name_angarms<<"hangle_arms"<<evtno;
			hangle_arms[good_evt] = new TH2F(name_angarms.str().c_str(),name_angarms.str().c_str(),45,0,180,45,-180,180);
			name_super<<"hsuperpose"<<evtno;
			hsuper[good_evt] = new TH2F(name_super.str().c_str(),name_angarms.str().c_str(),140,0,140,128,0,128);
			name_waveform<<"hwaveform"<<evtno;
			hwaveform[good_evt] = new TH2F(name_waveform.str().c_str(),name_waveform.str().c_str(),1200,0,1200,timetemp,0,timetemp);
			name_waveform<<"dc";
			hwaveform_dc[good_evt] = new TH2F(name_waveform.str().c_str(),name_waveform.str().c_str(),1200,0,1200,timetemp,0,timetemp);
			hwaveform_hit[good_evt] = new TH2F(Form("hwaveform_hit%d",evtno),Form("hwaveform_hit%d",evtno),1200,0,1200,timetemp,0,timetemp);
			name_energy<<"henergy"<<evtno;
			henergy[good_evt] = new TH1F(name_energy.str().c_str(),name_energy.str().c_str(),timetemp,0,timetemp);
			name_delta<<"hdeltat"<<evtno;
			hdeltat[good_evt] = new TH2F(name_delta.str().c_str(),name_delta.str().c_str(),140,0,140,128,0,128);
			name_flat<<"hflat"<<evtno;
			hflat[good_evt] = new TH2F(name_flat.str().c_str(),name_flat.str().c_str(),300,-150,150,300,-150,150);
			name_nn<<"hnearestneighbour"<<evtno;
			hnearestneighbour[good_evt] = new TH1F(name_nn.str().c_str(),name_nn.str().c_str(),50,0,50);
			name_arm1_xy<<"harm1_xy"<<evtno;
			harm1_xy[good_evt] = new TH2F(name_arm1_xy.str().c_str(),name_arm1_xy.str().c_str(),140,0,140,128,0,128);
			name_arm1_yz<<"harm1_yz"<<evtno;
			harm1_yz[good_evt] = new TH2F(name_arm1_yz.str().c_str(),name_arm1_yz.str().c_str(),140,0,140,timetemp,0,timetemp);
			name_mult<<"hmults"<<evtno;
			hmults[good_evt] = new TH1F(name_mult.str().c_str(),name_mult.str().c_str(),timetemp,0,timetemp);
			name_HoughXY<<"hHoughXY"<<evtno;
			hHoughXY[good_evt] = new TH2F(name_HoughXY.str().c_str(),name_HoughXY.str().c_str(),1800,0,180,600,-50,50);
			name_HoughXZ<<"hHoughXZ"<<evtno;
			hHoughXZ[good_evt] = new TH2F(name_HoughXZ.str().c_str(),name_HoughXZ.str().c_str(),1800,0,180,600,-50,50);
			name_HoughYZ<<"hHoughYZ"<<evtno;
			hHoughYZ[good_evt] = new TH2F(name_HoughYZ.str().c_str(),name_HoughYZ.str().c_str(),1800,0,180,600,-50,50);
                        name_HoughXY2<<"hHoughXY2_"<<evtno;
                        hHoughXY2[good_evt] = new TH2F(name_HoughXY2.str().c_str(),name_HoughXY2.str().c_str(),1800,0,180,600,-200,200);
                        name_HoughXZ2<<"hHoughXZ2_"<<evtno;
                        hHoughXZ2[good_evt] = new TH2F(name_HoughXZ2.str().c_str(),name_HoughXZ2.str().c_str(),1800,0,180,600,-200,200);
                        name_HoughYZ2<<"hHoughYZ2_"<<evtno;
//                        hHoughYZ2[good_evt] = new TH2F(name_HoughYZ2.str().c_str(),name_HoughYZ2.str().c_str(),1800,0,180,600,-200,200);
			name_distxy<<"hdistxy"<<evtno;
			hdistxy[good_evt] = new TH2F(name_distxy.str().c_str(),name_distxy.str().c_str(),140,0,245,128,0,224);
//			hdistxy[good_evt] = new TH2F(name_distxy.str().c_str(),name_distxy.str().c_str(),140,-100,100,128,-100,100);
			name_distyz<<"hdistyz"<<evtno;
			hdistyz[good_evt] = new TH2F(name_distyz.str().c_str(),name_distyz.str().c_str(),128,0,224,timetemp,0,256);
//			hdistyz[good_evt] = new TH2F(name_distyz.str().c_str(),name_distyz.str().c_str(),128,-100.,100.,timetemp,-100.,100.);
			name_distxz<<"hdistxz"<<evtno;
			hdistxz[good_evt] = new TH2F(name_distxz.str().c_str(),name_distxz.str().c_str(),140,0,245,timetemp,0,128);
//			hdistxz[good_evt] = new TH2F(name_distxz.str().c_str(),name_distxz.str().c_str(),140,-100,100,timetemp,-100,100.);

			hdistxy1[good_evt] = new TH2F(Form("hdxy1_%d",evtno),"",140,0,245,128,0,224);
			hdistxy2[good_evt] = new TH2F(Form("hdxy2_%d",evtno),"",140,0,245,128,0,224);
			hdistxy2a[good_evt] = new TH2F(Form("hdxy2a_%d",evtno),"",140,0,245,128,0,224);
			hdistxy2b[good_evt] = new TH2F(Form("hdxy2b_%d",evtno),"",140,0,245,128,0,224);
			hdistxy2c[good_evt] = new TH2F(Form("hdxy2c_%d",evtno),"",140,0,245,128,0,224);
			hdistxy3[good_evt] = new TH2F(Form("hdxy3_%d",evtno),"",140,0,245,128,0,224);
			hdistyz1[good_evt] = new TH2F(Form("hdyz1_%d",evtno),"",128,0,224,timetemp,0,256);
			hdistyz2[good_evt] = new TH2F(Form("hdyz2_%d",evtno),"",128,0,224,timetemp,0,256);
			hdistyz2a[good_evt] = new TH2F(Form("hdyz2a_%d",evtno),"",128,0,224,256,0,256);
			hdistyz2b[good_evt] = new TH2F(Form("hdyz2b_%d",evtno),"",128,0,224,256,0,256);
			hdistyz2c[good_evt] = new TH2F(Form("hdyz2c_%d",evtno),"",128,0,224,256,0,256);
			hdistyz3[good_evt] = new TH2F(Form("hdyz3_%d",evtno),"",128,0,224,256,0,256);
			hdistxz1[good_evt] = new TH2F(Form("hdxz1_%d",evtno),"",140,0,245,256,0,256);
			hdistxz2[good_evt] = new TH2F(Form("hdxz2_%d",evtno),"",140,0,245,256,0,256);
			hdistxz2a[good_evt] = new TH2F(Form("hdxz2a_%d",evtno),"",140,0,245,256,0,256);
			hdistxz2b[good_evt] = new TH2F(Form("hdxz2b_%d",evtno),"",140,0,245,256,0,256);
			hdistxz2c[good_evt] = new TH2F(Form("hdxz2c_%d",evtno),"",140,0,245,256,0,256);
			hdistxz3[good_evt] = new TH2F(Form("hdxz3_%d",evtno),"",140,0,245,256,0,256);
			harmxy[good_evt] = new TH2F(Form("harmxy_%d",evtno),"",140,0,245,128,0,224);
			harmyz[good_evt] = new TH2F(Form("harmyz_%d",evtno),"",128,0,224,512,0,256);
			harmxz[good_evt] = new TH2F(Form("harmxz_%d",evtno),"",140,0,245,512,0,256);
			name_projXY<<"hHoughprojXY"<<evtno;
			hHoughprojXY[good_evt] = new TH1F(name_projXY.str().c_str(),name_projXY.str().c_str(),1800,0,1800);
			name_projYZ<<"hHoughprojYZ"<<evtno;
			hHoughprojYZ[good_evt] = new TH1F(name_projYZ.str().c_str(),name_projYZ.str().c_str(),1800,0,1800);
			name_projXZ<<"hHoughprojXZ"<<evtno;
			hHoughprojXZ[good_evt] = new TH1F(name_projXZ.str().c_str(),name_projXZ.str().c_str(),1800,0,1800);
			name_stage1<<"hstage1_"<<evtno;
			hstage1[good_evt] = new TH2F(name_stage1.str().c_str(),name_stage1.str().c_str(),140,-50,50,128,-50,50);
			name_stage2<<"hstage2_"<<evtno;
			hstage2[good_evt] = new TH2F(name_stage2.str().c_str(),name_stage2.str().c_str(),140,-50,50,128,-50,50);
			name_stage3<<"hstage3_"<<evtno;
			hstage3[good_evt] = new TH2F(name_stage3.str().c_str(),name_stage3.str().c_str(),140,-50,50,128,-50,50);
			name_chains<<"hchainsevt"<<evtno;
			hchainsevt[good_evt] = new TH1F(name_chains.str().c_str(),name_chains.str().c_str(),timetemp,0,timetemp);
			name_maxchains<<"hmaxchains"<<evtno;
                        hmaxchains[good_evt] = new TH2F(name_maxchains.str().c_str(),name_maxchains.str().c_str(),140,0,140,128,0,128);
                        name_maxstrips<<"hmaxstrips"<<evtno;
                        hmaxstrips[good_evt] = new TH2F(name_maxstrips.str().c_str(),name_maxstrips.str().c_str(),140,0,140,128,0,128);

		}
	return;
}


double Analysis::RemoveBG() {
	Float_t waveformBG[512];//array to hold the background for the waveform
	int check_all_sat=0;
	integratewf=0.;
	if(mmCobo[waveformno]==1 && mmAsad[waveformno]==0 && mmAget[waveformno]==0) {//select si fronts
		for(int t=0;t<timetemp;t++) mmWaveformY[waveformno][t]=4095-mmWaveformY[waveformno][t];//flip the silicon fronts
	}
	double backgrnd = 1000;
	TH1F *hfindbg = new TH1F("hfindbg","Find BG via mode of waveform",2048,0,2048);
	saturated=false;
	double average_start=0,average_end=0;
	int start_no=0,end_no=0;
	for(int t=0;t<timetemp;t++) {
		hfindbg->Fill(mmWaveformY[waveformno][t]);
		if(mmWaveformY[waveformno][t]==4095) {
			check_all_sat++;//this timebucket is saturated
			saturated=true;
		}
		if(t<40) {
			average_start+=mmWaveformY[waveformno][t];
			start_no++;
		}
		if(t>511-40) {
			average_end+=mmWaveformY[waveformno][t];
			end_no++;
		}
	}
	average_start/=1.*start_no;
	average_end/=1.*end_no;
	if(average_start-average_end>80) {
		step_count++;
	}
	double backgrnd_E=hfindbg->GetMaximumBin()-baseline_set;//get mode
//	double backgrnd_E=hfindbg->GetMaximumBin()-50;//get mode
//	backgrnd_E=0.;
//	cout<<baseline_set<<endl;
	if(check_all_sat>100) backgrnd_E=4095;
	delete hfindbg;
	for(int t=0;t<timetemp;t++) {//integrate up until the last peaking waveform
		waveformBG[t]=mmWaveformY[waveformno][t];
//		totalEnergy+=mmWaveformY[waveformno][t]-backgrnd_E-baseline_set;//count total energy as integral of waveform
		integratewf+=mmWaveformY[waveformno][t]-backgrnd_E-baseline_set;
		hbg->SetBinContent(t,mmWaveformY[waveformno][t]);
		if(mmWaveformY[waveformno][t]<backgrnd) backgrnd=mmWaveformY[waveformno][t];//take minimum, normal histo method doesn't seem to work very well
		wave[t]=double(mmWaveformY[waveformno][(t-source_offset+timetemp)%timetemp]);

	}
//	cout<<backgrnd_E<<endl;
	hwave = new TH1F("hwave","Current waveform",timetemp,0,timetemp);
	hwave2 = new TH1F("hwave2","Current waveform",timetemp,0,timetemp);
	for(int t=0;t<10;t++) mmWaveformY[waveformno][t]=backgrnd_E;//get rid of the bad start of the waveform that affects the convolution (set to backgrnd_E so when removed ==0)
	for(int t=500;t<512;t++) mmWaveformY[waveformno][t]=backgrnd_E;//get rid of the bad start of the waveform that affects the convolution (set to backgrnd_E so when removed ==0)
	double maxE=0;
	for(int t=0;t<timetemp;t++) {//CORRECTED WAVEFORM
		if(good_evt<draw_events) {
//			hpulse[good_evt]->Fill(t,mmWaveformY[waveformno][t]-backgrnd_E,map_y[mmAsad[waveformno]][mmAget[waveformno]][mmChan[waveformno]]+2);//draw all waveforms - colour determined by distance along the pad
//			cout<<good_evt<<"\t"<<t<<"\t"<<mmWaveformY[waveformno][t]<<"\t"<<backgrnd_E<<map_y[mmAsad[waveformno]][mmAget[waveformno]][mmChan[waveformno]]+2<<endl;
			henergy[good_evt]->Fill(t,mmWaveformY[waveformno][t]-backgrnd_E-baseline_set);
//			hpulse[good_evt]->Fill(mmWaveformY[waveformno][t],40);

		}
//		cout<<t<<"\t"<<mmWaveformY[waveformno][t];
		mmWaveformY[waveformno][t]-=backgrnd_E;//take off the background
		wave[t]=mmWaveformY[waveformno][t];
//		cout<<t<<"\t"<<mmWaveformY[waveformno][t]<<endl;
//		if(t>source_offset)hwave->SetBinContent(t,mmWaveformY[waveformno][(t-source_offset+timetemp)%timetemp]);
//		else hwave->SetBinContent(t,0);
		hwave->SetBinContent(t,mmWaveformY[waveformno][(t-source_offset+timetemp)%timetemp]);
		hwave2->SetBinContent(t,mmWaveformY[waveformno][t]);
//		wave[t]=double(mmWaveformY[waveformno][t]);
//		cout<<t<<"\t"<<wave[t]<<"\t"<<source_offset<<"\t"<<(t-source_offset+timetemp)%timetemp<<"\t"<<mmWaveformY[waveformno][t]<<endl;
		if(good_evt<draw_events) hwaveform[good_evt]->Fill(waveformno,t,wave[t]);
//		cout<<good_evt<<"\t"<<waveformno<<"\t"<<t<<"\t"<<wave[t]<<endl;
//		if(t<=source_offset) wave[t]=0.;
		if(wave[t]>maxE) maxE=wave[t];

	}
	maxE-=baseline_set;
	if(average_start-average_end>80) maxE=0;
	return maxE;
}


double Analysis::GetEnergy4He(double range) {
//	return -8.597e-5*range*range+4.2556e-2*range-2.3766e-1;
	if(pressure==100) range*=2.;//Convert from the 100 Torr range to the 50 Torr range
	double end=0.1;
	return -2.5156e-5*range*range+1.6060e-2*range-1.2113e-1+end;//New range based on GEANT4

//	return -1.8045e-5*range*range+2.0237e-2*range-1.8144e-1;
}

double Analysis::GetEnergy13C(double range) {//50 Torr
	if(pressure==100) range*=2.;

	return 0.0052*TMath::Power(range,1.5641);
}

double Analysis::GetEnergy9Be(double range) {
	if(pressure==100) range*=2.;

	return 0.0034*TMath::Power(range,1.5462);//MeV
}

double Analysis::GetEnergy12C(double range) {
	if(pressure==100) range*=2.;

	return 0.0053*TMath::Power(range,1.5677);//MeV from mm
}

void Analysis::AnalyzeSi() {
	int m=waveformno;
	if(mmCobo[m]==0) return;
	int select=0;
	if(mmChan[m]>=34 && mmChan[m]<51 && mmAsad[m]==0) select=1;
	if(mmChan[m]>=34 && mmChan[m]<38 && mmAsad[m]==1) select=2;
	if(select==0) return;
//	Only selecting relevant front back channels to sweep through
	double maxE=0,maxT=0;
	for(int t=0;t<timetemp;t++) {
		if(mmWaveformY[m][t]>maxE) {
			maxE=mmWaveformY[m][t];//new max
			maxT=t;
		}
	}
	int silicon_channel=0;
	silicon_channel=mmChan[m]-34+mmAsad[m]*16;//0->15 are fronts then 16,17,18,19 are backs
	hSilicon[silicon_channel]->Fill(maxE);
	return;
}


void Analysis::ReadCalibrationFile(TString cname) {

	ifstream incalib;
	incalib.open(cname);
	if(!incalib) cout<<BOLDRED<<"Thickness file not found\n!"<<RESET;
	if(incalib) cout<<BOLDGREEN<<"Thicknesses read\n"<<RESET;
	int ccount=0;
	double thisgain,thisoffset;
	string channel;
	while(incalib.good() && ccount<20){
        	incalib>>thisgain>>thisoffset;
	        gain[ccount]=thisgain;
	        offset[ccount]=thisoffset;;
		ccount++;
//		cout<<ccount<<endl;
	}

	return;
}



double Analysis::CalcChisqrdDist(const double *par) {
//      Chi-squared (or thing to be minimized) is the average distance from each voxel hit to the NEAREST line of the 3 alphas (weighted by enery deposited in the voxel). The free parameters are the end$
        Double_t chisqrd=0,totaldist=0.,accumulatedE=0.;
	double bestdiff=10000;
	double arm_length[2]={0.};
	double distance_[2]={1000};
	bool debug_chi=false;
	if(server_mode) debug_chi=false;
	Eright=0;
	Eleft=0.;
	double maxheight=0;
	double maxy=0,miny=224;
	double maxx=0,minx=224;
	double maxz=0,minz=512;
	double dist_end[2]={100000.,100000.};
	double dist_vertex=10000.;
	double furthest_p[2]={0.};//total point length of each arm
	int furthest_point[2]={0};//total point index of each arm
	double Eweight=0;
        for(unsigned int numi=0;numi<pointx.size();numi++) {
                double tmpx,tmpy,tmpz;
                tmpx=pointx.at(numi);
                tmpy=pointy.at(numi);
                tmpz=pointz.at(numi);
		Eweight+=(pointE.at(numi)*pointE.at(numi));
                double x0[3]={tmpx,tmpy,tmpz};//current point
                double x1[3]={par[0],par[1],par[2]};//center of decay
//                                              Distance between line between x1 and x2 and the point x0 is:
//                                              d=|(x2-x1)X(x1-x0)|/|x2-x1|
		if(tmpy<miny) miny=tmpy;
		if(tmpy>maxy) maxy=tmpy;
		if(tmpx<minx) minx=tmpx;
		if(tmpx>maxx) maxx=tmpx;
		if(tmpz>maxz) maxz=tmpz;
		if(tmpz<minz) minz=tmpz;
		if(maxz-minz>maxheight) maxheight=maxz-minz;

                double cross1x,cross1y,cross1z,bestdiff=1000.;
		bestdiff=10000;
		if(x0[0]<x1[0]) Eleft+=pointE.at(numi);
		if(x0[0]>x1[0]) Eright+=pointE.at(numi);
                for(int arm=0;arm<2;arm++) {//do for each arm
                        double x2[3]={par[arm*3+3],par[arm*3+4],par[arm*3+5]};
//                                              cout<<x2[0]<<"\t"<<x2[1]<<"\t"<<x2[2]<<endl;
                        Cross(x2[0]-x1[0],x2[1]-x1[1],x2[2]-x1[2],x1[0]-x0[0],x1[1]-x0[1],x1[2]-x0[2],cross1x,cross1y,cross1z);
                        double num=sqrt(cross1x*cross1x+cross1y*cross1y+cross1z*cross1z);
                        double denom=sqrt(TMath::Power(x2[0]-x1[0],2)+TMath::Power(x2[1]-x1[1],2)+TMath::Power(x2[2]-x1[2],2));
                        double thisdist=num/denom;
                        double distcentr=sqrt(TMath::Power(x1[0]-x0[0],2)+TMath::Power(x1[1]-x0[1],2)+TMath::Power(x1[2]-x0[2],2));//distance from this point to interaction
                        double disttoend=sqrt(TMath::Power(x2[0]-x0[0],2)+TMath::Power(x2[1]-x0[1],2)+TMath::Power(x2[2]-x0[2],2));//distance from the end to point
                        double distarm=sqrt(TMath::Power(x2[0]-x1[0],2)+TMath::Power(x2[1]-x1[1],2)+TMath::Power(x2[2]-x1[2],2));//distance from the end to the interaction
			if(disttoend<dist_end[arm]) dist_end[arm]=disttoend;
                        double d_point[3],d_arm[3],dpointtoarm[3];
                        double dist_to_centr=0.,dist_to_arm=0.,dist_ptoarm=0.;
			distance_[arm]=thisdist;
			if(distcentr>furthest_p[arm] && thisdist<5) {
				furthest_p[arm]=distcentr;//Find the furthest points for each arm
				furthest_point[arm]=numi;
			}
			if(distcentr>distarm) thisdist=10.*disttoend;//If the point is PAST the arm then instead of the perpendicular distance we instead have to use the distance from the point to the end of the arm with a 10x penalty
			if(thisdist<bestdiff) {
				bestdiff=thisdist;
			}
		}
		double distance_to_vertex=sqrt(TMath::Power(x1[0]-x0[0],2)+TMath::Power(x1[1]-x0[1],2)+TMath::Power(x1[2]-x0[2],2));
		if(distance_to_vertex<dist_vertex) dist_vertex=distance_to_vertex;
		int closest=0;
		if(distance_[0]>distance_[1]) {
			closest=1;
		}
		if(distance_to_vertex>arm_length[closest]) {
			arm_length[closest]=distance_to_vertex;
		}
		chisqrd+=(bestdiff*bestdiff*pointE.at(numi)*pointE.at(numi))/(1.75*1.75);//position error of 1.75 mm assumed
//		cout<<bestdiff<<"\t"<<chisqrd<<"\t"<<chisqrd/(num-9)<<"\t"<<num<<"\t\t"<<tmpx<<"\t"<<tmpy<<"\t"<<tmpz<<endl;
	}
	chisqrd/=(Eweight*(double(pointx.size())-9.));//NDOF for normal fitting
	if(debug_chi) cout<<"Fitting chisqrd = "<<chisqrd;
	double xzangle1=atan2(par[5],par[3]);
	double xzangle2=atan2(par[8],par[6]);
	if(good_evt<draw_events) hdistyz2a[good_evt]->SetTitle(Form("Left/right E    %g %g",Eleft,Eright));
	chisqrd+=(0.01*dist_end[0]*dist_end[0]+0.01*dist_end[1]*dist_end[1]+0.2*dist_vertex*dist_vertex);//Ensure the fitted ends align with the data
	double dparm[2]={0};//distance from end of two arms to the fitted ends
//	cout<<furthest_point[0]<<"\t"<<furthest_point[1]<<"\t"<<pointx.size()<<endl;
	for(int arm=0;arm<2;arm++) {
		dparm[arm]=sqrt(TMath::Power(pointx.at(furthest_point[arm])-par[arm*3+3],2)+TMath::Power(pointy.at(furthest_point[arm])-par[arm*3+4],2)+TMath::Power(pointz.at(furthest_point[arm])-par[arm*3+5],2));
	}
//	chisqrd+=(dparm[0]*dparm[0]+dparm[1]*dparm[1])*0.1;
	if(debug_chi) cout<<"Arm contributions = "<<(dparm[0])<<"\t"<<(dparm[1])<<endl;
	if(debug_chi)cout<<"End contributions = "<<(dist_end[0])/1.75<<"\t,\t"<<(dist_end[1])*1.75<<"\t,\t"<<dist_vertex/1.75<<endl;
	if(debug_chi) cout<<"Total track chisqrd = "<<chisqrd<<endl;
//	chisqrd+=(fabs(xzangle1-xzangle2-4.*atan(1.)))/0.1;
//	if(debug_chi) cout<<"Angle penalty: "<<(fabs(xzangle1-xzangle2-4.*atan(1.)))/0.1<<endl;
	if(par[3]-par[0]>0 && par[6]-par[0]>0) {
		chisqrd+=10000*(par[3]-par[0])*(par[6]-par[0]);//make opposite sides in x
		if(debug_chi) cout<<"Opposite x sides penalty"<<endl;
	}
	if(par[3]-par[0]<0 && par[6]-par[0]<0) {
		chisqrd+=10000*(par[0]-par[3])*(par[0]-par[6]);//make opposite sides
		if(debug_chi) cout<<"Opposite x sides penalty"<<endl;
	}
	if(par[5]-par[2]>0 && par[8]-par[2]>0) {
		chisqrd+=10000*(par[5]-par[2])*(par[8]-par[2]);//make opposite sides in z
		if(debug_chi) cout<<"Opposite z sides penalty"<<endl;
	}
	if(par[5]-par[2]<0 && par[8]-par[2]<0) {
		chisqrd+=10000*(par[2]-par[5])*(par[2]-par[8]);//make opposite sides
		if(debug_chi) cout<<"Opposite z sides penalty"<<endl;
	}
	if(fabs(par[0]-122.5)>15) {
		chisqrd+=(fabs(par[0]-122.5)-15)*(fabs(par[0]-122.5)-15)*1e5;//make it vaguely central
		if(debug_chi) cout<<"Non central penalty"<<endl;
	}
//	Calculate how good the angles are from theta12C,psi12C
	double smallest_th_dist[2]={180,180.};
	double totaldx1=-par[0]+par[3];
	double totaldx2=-par[0]+par[6];
	double totaldy1=-par[1]+par[4];
	double totaldy2=-par[1]+par[7];
	double totaldz1=-par[2]+par[5];
	double totaldz2=-par[2]+par[8];
	double theta_here=atan2(sqrt(totaldx1*totaldx1+totaldz1*totaldz1),totaldy1)/d2r;
	double psi_here=atan2(sqrt(totaldx2*totaldx2+totaldz2*totaldz2),totaldy2)/d2r;
	int heavyis=0;
	double heavyE=0,lightE=0;
	double energy_contr[2];
	TString debugstring[2];
//	double ftheta[2][2]={0.};//for two light/heavy combos for theta/psi
//	double dE_heavy=0.,dE_light=0.;
	double dE_track_heavy[3]={0.};
	double dE_track_light[3]={0.};
//	for(int looper=0;looper<2;looper++) {//loops over interchanging light and heavy
	debugstring[looper]="debug looper ";
	debugstring[looper]+=looper;
	double diste1=dist_end[0];
	double diste2=dist_end[1];
	debugstring[looper]+=" distends = ";
	debugstring[looper]+=diste1;
	debugstring[looper]+=" , ";
	debugstring[looper]+=diste2;
	debugstring[looper]+=" , ";
	debugstring[looper]+=dist_vertex;
	for(int i=0;i<theta12C.size();i++) {
		double dist=1e9;
		if(looper==0) {//For looper 0 assume that theta [LEFT] is matched to theta theory
			if(select_nalpha==1)dist=sqrt(TMath::Power(theta_here-theta12C.at(i),2)+TMath::Power(psi_here-psi12C.at(i),2));
			if(select_nalpha==2)dist=sqrt(TMath::Power(theta_here-theta16O.at(i),2)+TMath::Power(psi_here-psi16O.at(i),2));
			if(select_nalpha==3)dist=sqrt(TMath::Power(theta_here-theta16Op.at(i),2)+TMath::Power(psi_here-psi16Op.at(i),2));
			if(dist<smallest_th_dist[looper]) {
				smallest_th_dist[looper]=dist;
				fit_theta[select_nalpha-1][0]=theta_here;
				fit_theta[select_nalpha-1][1]=psi_here;
				ftheta[looper][0]=theta_here;
				ftheta[looper][1]=psi_here;
				if(select_nalpha==1) {
					theory_theta[0]=theta12C.at(i);
					theory_theta[1]=psi12C.at(i);
					heavyE=E12C_heavy.at(i);
					lightE=E12C_light.at(i);
				}
				if(select_nalpha==2) {
					theory_theta[0]=theta16O.at(i);
					theory_theta[1]=psi16O.at(i);
					heavyE=E16O_heavy.at(i);
					lightE=E16O_light.at(i);
				}
				if(select_nalpha==3) {
					theory_theta[0]=theta16Op.at(i);
					theory_theta[1]=psi16Op.at(i);
					heavyE=E16Op_heavy.at(i);
					lightE=E16Op_light.at(i);
				}
				heavyis=1;
//			Set global variables here for angle?
			}
		}
		if(looper==1) {//For looper 1 assume that theta is matched to psi
			if(select_nalpha==1)dist=sqrt(TMath::Power(theta_here-psi12C.at(i),2)+TMath::Power(psi_here-theta12C.at(i),2));//test other way around
			if(select_nalpha==2)dist=sqrt(TMath::Power(theta_here-psi16O.at(i),2)+TMath::Power(psi_here-theta16O.at(i),2));//test other way around
			if(select_nalpha==3)dist=sqrt(TMath::Power(theta_here-psi16Op.at(i),2)+TMath::Power(psi_here-theta16Op.at(i),2));//test other way around
			if(dist<smallest_th_dist[looper]) {
				smallest_th_dist[looper]=dist;
//		Set global variables here for angle?
				fit_theta[select_nalpha-1][0]=theta_here;
				fit_theta[select_nalpha-1][1]=psi_here;
				ftheta[looper][0]=psi_here;
				ftheta[looper][1]=theta_here;
					if(select_nalpha==1) {
					theory_theta[1]=theta12C.at(i);
					theory_theta[0]=psi12C.at(i);
					heavyE=E12C_heavy.at(i);
					lightE=E12C_light.at(i);
				}
				if(select_nalpha==2) {
					theory_theta[1]=theta16O.at(i);
					theory_theta[0]=psi16O.at(i);
					heavyE=E16O_heavy.at(i);
					lightE=E16O_light.at(i);
				}
				if(select_nalpha==3) {
					theory_theta[1]=theta16Op.at(i);
					theory_theta[0]=psi16Op.at(i);
					heavyE=E16Op_heavy.at(i);
					lightE=E16Op_light.at(i);
				}
				heavyis=0;
			}
		}
	}
//	Get ranges
	int lightis=0;
	if(heavyis==0) lightis=1;
	double dx_heavy=0,dy_heavy=0,dz_heavy=0;
	double dx_light=0,dy_light=0,dz_light=0;
	dx_heavy=par[3*heavyis+3]-par[0];
	dy_heavy=par[3*heavyis+4]-par[1];
	dz_heavy=par[3*heavyis+5]-par[2];
//		dx_light=par[3*lightis+3]-par[0];
	dy_light=par[3*lightis+4]-par[1];
	dz_light=par[3*lightis+5]-par[2];

	bool light_escapes=false,heavy_escapes=false;
	if(fabs(dz_light)>50 && maxheight>50) {//If data looks like it escapes and the fit represents the line going past this then escape
		light_escapes=true;
	}
	if(fabs(dz_heavy)>50 && maxheight>50) {
		heavy_escapes=true;
	}
	if((par[3*lightis+4]>210 && maxy>200) || (miny<10 && par[3*lightis+4]<0.)) {//escapes off the start/end of MM
		light_escapes=true;
	}
	if((par[3*heavyis+4]>210 && maxy>200) || (miny<10 && par[3*heavyis+4]<0.)) {//escapes off the start/end of MM
		heavy_escapes=true;
	}
	if((par[3*lightis+3]>224 && maxx>210) || (minx<10 && par[3*lightis+3]<0.)) {//escapes off the side of MM
		light_escapes=true;
	}
	if((par[3*heavyis+3]>224 && maxx>210) || (minx<10 && par[3*heavyis+3]<0.)) {//escapes off the side of MM
		heavy_escapes=true;
	}
	if(light_escapes) debugstring[looper]+=" light escapes ";
	if(heavy_escapes) debugstring[looper]+=" heavy escapes ";
	double R_heavy,R_light;
	R_heavy=sqrt(dx_heavy*dx_heavy+dy_heavy*dy_heavy+dz_heavy*dz_heavy);
	R_light=sqrt(dx_light*dx_light+dy_light*dy_light+dz_light*dz_light);
	dE_heavy=0.;
	dE_light=0.;

	if(select_nalpha==1) {
		dE_heavy=fabs(GetEnergy9Be(R_heavy)-heavyE);//12C
//		cout<<"Heavy 12C: "<<GetEnergy9Be(R_heavy)<<"\t"<<heavyE<<endl;
		if(heavy_escapes /*&& heavyE<GetEnergy9Be(R_heavy)*/) dE_heavy=0.02;//small penalty for escaping
//		if(!heavy_escapes)hheavyR->Fill(GetEnergy9Be(R_heavy),heavyE);
		if(!heavy_escapes && GetEnergy9Be(R_heavy)-heavyE) {
			rheavy[select_nalpha-1]=GetEnergy9Be(R_heavy);
			rtheavy[select_nalpha-1]=heavyE;
		}

	}
	if(select_nalpha==2) {
		dE_heavy=fabs(GetEnergy13C(R_heavy)-heavyE);//16O
//		cout<<"Heavy 16O: "<<GetEnergy9Be(R_heavy)<<"\t"<<heavyE<<endl;
		if(heavy_escapes /*&& heavyE<GetEnergy13C(R_heavy)*/) dE_heavy=0.02;//small penalty for escaping
//		if(!heavy_escapes)hheavyR->Fill(GetEnergy13C(R_heavy),heavyE);
		if(!heavy_escapes && GetEnergy13C(R_heavy)<heavyE) {
			rheavy[select_nalpha-1]=GetEnergy13C(R_heavy);
			rtheavy[select_nalpha-1]=heavyE;
		}

	}
	if(select_nalpha==3) {
		dE_heavy=fabs(GetEnergy13C(R_heavy)-heavyE);//16Op
//		cout<<"Heavy 16O: "<<GetEnergy9Be(R_heavy)<<"\t"<<heavyE<<endl;
		if(heavy_escapes /*&& heavyE<GetEnergy13C(R_heavy)*/) dE_heavy=0.02;//small penalty for escaping
//		if(!heavy_escapes)hheavyR->Fill(GetEnergy13C(R_heavy),heavyE);
		if(!heavy_escapes && GetEnergy13C(R_heavy)<heavyE) {
			rheavy[select_nalpha-1]=GetEnergy13C(R_heavy);
			rtheavy[select_nalpha-1]=heavyE;
		}

	}
//	if(!light_escapes)hlightR->Fill(GetEnergy4He(R_light),lightE);
	if(!light_escapes) {
		rlight[select_nalpha-1]=GetEnergy4He(R_light);
		rtlight[select_nalpha-1]=lightE;
	}
	dE_light=fabs(GetEnergy4He(R_light)-lightE);//alphas for both
	if(light_escapes /*&& lightE<GetEnergy4He(R_light)*/) dE_light=0.02;//No huge penalty for ends not meeting up -- modest one to avoid simple moving out of bounds
//	cout<<"Light: "<<GetEnergy4He(R_light)<<"\t"<<lightE<<endl;
	if(debug_chi) {
		cout<<"Light escape: "<<light_escapes<<"\tHeavy escapes: "<<heavy_escapes<<endl;
		cout<<"dE_light: "<<dE_light<<"\tdE_heavy: "<<dE_heavy<<endl;
	}
	energy_contr[looper]=(dE_light*dE_light+dE_heavy*dE_heavy);
	debugstring[looper]+=" dE_light = ";
	debugstring[looper]+=(dE_light);
	debugstring[looper]+=" dE_heavy = ";
	debugstring[looper]+=(dE_heavy);
	dE_track_light[looper]=dE_light;
	dE_track_heavy[looper]=dE_heavy;
//	}
	if(good_evt<draw_events) hdistyz2b[good_evt]->SetTitle(debugstring[0]);
	if(good_evt<draw_events) hdistxz2b[good_evt]->SetTitle(debugstring[1]);
	int sellooper=looper;
	dlength[select_nalpha-1][0]=dE_track_light[sellooper];
	dlength[select_nalpha-1][1]=dE_track_heavy[sellooper];
	dtheta[select_nalpha-1]=smallest_th_dist[sellooper];
	if(debug_chi) cout<<"Current chisqrd = "<<chisqrd<<endl;
	chisqrd+=(smallest_th_dist[sellooper]*smallest_th_dist[sellooper])/1000;//Make a loose chi-sqrd optimization for angle dtheta=3
	if(debug_chi) cout<<"Angle contr: "<<(smallest_th_dist[sellooper]*smallest_th_dist[sellooper])/4.<<endl;
//	chisqrd+=energy_contr[sellooper];//add the best heavy/light energy match - sigma of 0.1 MeV
//	chisqrd+=dE_heavy*dE_heavy;
	if(debug_chi) cout<<"Energy contr: "<<50*energy_contr[sellooper]<<endl;
	fit_theta[select_nalpha-1][0]=ftheta[sellooper][0];
	fit_theta[select_nalpha-1][1]=ftheta[sellooper][1];
	if(debug_chi) {
		for(int i=0;i<9;i++) {
			cout<<"Par "<<i<<"\t"<<par[i]<<endl;
		}
	}
//	chisqrd/=5;//turn back to reduced chi-sqrd
	if(debug_chi)cout<<"Nalpha chi= "<<chisqrd<<" for "<<select_nalpha<<endl;
//	chisqrd_nalpha=(smallest_th_dist[sellooper]*smallest_th_dist[sellooper])/25.+(energy_contr[sellooper])*20;
	chisqrd_nalpha_theta=(smallest_th_dist[sellooper]*smallest_th_dist[sellooper]);
	chisqrd_nalpha=(smallest_th_dist[sellooper]*smallest_th_dist[sellooper])/16.+dE_heavy*dE_heavy/(0.4*0.4);
	if(pressure==100) chisqrd_nalpha=(smallest_th_dist[sellooper]*smallest_th_dist[sellooper])/64.+dE_heavy*dE_heavy/(0.5*0.5);
	chisqrd_nalpha_dist=(energy_contr[sellooper]);
//	chisqrd_nalpha=chisqrd;
	if(debug_chi) cout<<"Evaluated chisqrd = "<<chisqrd_nalpha<<endl;
//	cout<<select_nalpha<<"\t"<<looper<<"\t"<<chisqrd<<endl;
	dE_light_[looper]=dE_light;
	dE_heavy_[looper]=dE_heavy;
	return chisqrd;
}
