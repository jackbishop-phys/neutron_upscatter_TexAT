//8.25 MeV roughly centroid neutrons
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TGraphSmooth.h"
#include "TMath.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TH1.h"
#include <iostream>
#include <fstream>
void overlay() {
	ifstream in1;
	in1.open("gas_303/tag0030.xy");
//	in1.open("dsk3001/tag0097.xy");

	ifstream in2;
	in2.open("no_gas_304/tag0030.xy");
	ifstream in3;
	in3.open("gas_303/tag0031.xy");
//	in3.open("dsk3001/tag0098.xy");

	double time[4096],gamma[4096],neutron[4096],gammabg[4096],gamma2[4096];
	int count=0;
	while(in1.good()) {
		double timein,gammain;
		in1>>timein>>gammain;
		time[count]=timein;
		gamma[count]=gammain;
		count++;
	}
	count=0;
	while(in2.good()) {
		double timein,gammain;
		in2>>timein>>gammain;
		gammabg[count]=gammain;
		gamma2[count]=-gamma[count]+gammain;
		count++;	}
	count=0;
	while(in3.good()) {
		double timein,neutronin;
		in3>>timein>>neutronin;
//		if(count>150 && count<4096) {
		if(count>10 && count<8192) {
			neutron[count]=neutronin;
			cout<<count<<endl;
		}
		else {
			neutron[count]=0.;
		}
		count++;
	}
//
//	TGraph *g1 = new TGraph("gas_303/tag0030.xy");
//	TGraph *g2 = new TGraph("no_gas_304/tag0030.xy");
//	TGraph *g3 = new TGraph("gas_303/tag0031.xy");
	TGraph *g1 = new TGraph(count,time,gamma);
	TGraph *g2 = new TGraph(count,time,gammabg);
	TGraph *g3 = new TGraph(count,time,neutron);
	TGraph *g4 = new TGraph(count,time,gamma2);
	g1->SetLineColor(kBlack);
	g1->SetLineWidth(2);
	g1->SetMarkerColor(kBlack);
	g1->GetXaxis()->SetRange(100,300);
	g1->SetMaximum(300);
//	TCanvas *c1 = new TCanvas("c1","",0,0,800,600);
	TCanvas *c1 = new TCanvas();
	g1->Draw();
	g2->SetLineColor(kRed);
	g2->SetMarkerColor(kRed);
	g2->SetLineWidth(2);
//	g2->Draw("SAME");
	g4->SetLineColor(kBlue);
	g4->SetLineWidth(2);
	g4->SetMarkerColor(kBlue);
	g4->Draw("SAME");
	TLegend *leg = new TLegend(0.7,0.9,0.7,0.9);
	leg->AddEntry(g1,"Gas");
	leg->AddEntry(g2,"No gas");
	leg->Draw("SAME");
//	TCanvas *c2 = new TCanvas("c2","Neutron energy",0,600,800,600);
	TCanvas *c2 = new TCanvas();
	g3->Draw();
	g3->GetXaxis()->SetRangeUser(700,1000);
//	Convert to energies
//	E=0.5*m*v^2
//	v=d/t --> E=0.5*m*d^2/t^2
//	const int nbins=2000;
	const int nbins=1000;
	double Eneutron[nbins];
	double Eneutronnew[nbins];
	double errEneutronnew[nbins];
	double distance_havar = (29.904+0.0399+0.0254);//29.904 m @ + 4 cm to Havar+AVERAGE of 1" into NE213
	double distance_gamma = (29.904+0.034+0.0399+0.0254);//29.94 m @ speed of light + 2.4 cm offset for gas cell center + 4cm to Havar +  3.4 cm for 1/8" collimator (this is our gamma peak)+AVERAGE of 1 cm into NE213

	for(int i=0;i<nbins;i++) {
//		time[i]=((170-time[i]+3628)*0.262747)*1e-9;
//		Assume that the interaction occurs DIRECTLY after the Havar foil
		double time_gamma=distance_gamma/3e8;//ns for gamma ray to travel
//		Two gamma peaks are at 167.822 and 182.356
		cout<<time[i]<<endl;
		time[i]=((167.822-time[i])*0.216771)*1e-9+800e-9+time_gamma;//extra time + time diff between neutrons and gamma
//		time[i]=((1724-time[i])*0.216771)*1e-9+800e-9+time_gamma;//extra time + time diff between neutrons and gamma
//time[i]=((167-time[i]+(800./0.211))*0.211)*1e-9+time_gamma;//extra time + time diff between neutrons and gamma
//150 from gamma to beginning of TAC
		Eneutron[i]=0.5*1.67e-27*TMath::Power(distance_havar-0.0399,2.)/(time[i]*time[i]);//energy assuming interaction occurs at gas-cell centre
		Eneutron[i]/=1.6e-13;//J->MeV
		cout<<time[i]<<"\t"<<Eneutron[i]<<endl;
		double Edeuteron=1.0335*Eneutron[i]-3.5134;//From zero degrees scattering equation
//		Make the correction Carl suggested for time of flight for the deuteron
//		Work out WHERE in the gas-cell this interaction occured
		double distance_in_gas=0;
		double fraction=1.-(Edeuteron-5.)/0.21;//5.21->5 MeV for 0->1 fraction
		distance_in_gas=fraction*7.98*1e-2;//
//		cout<<Eneutron[i]<<"\t"<<Edeuteron<<"\t"<<fraction<<endl;
		double neutron_distance=distance_havar-distance_in_gas;
		double time_for_beam=distance_in_gas/sqrt((2.*Edeuteron*1.6e-13)/(1.67e-27));//t=d/v
//		cout<<"Time for beam "<<time_for_beam<<endl;
		double neutron_cell_time=distance_in_gas/sqrt((2.*Eneutron[i]*1.6e-13)/(1.67e-27));
		double time_for_neutron=time[i]-time_for_beam;
		Eneutronnew[i]=(0.5*1.67e-27*TMath::Power(neutron_distance,2.)/(time_for_neutron*time_for_neutron))/1.6e-13;
		errEneutronnew[i]=0.;
		if(i>920) Eneutronnew[i]=Eneutron[i];
//		cout<<i<<"\t"<<Eneutronnew[i]<<endl;
//		cout<<"energy shift "<<Eneutronnew[i]-Eneutron[i]<<"\t"<<fraction<<"\t"<<time_for_beam<<endl;
//		if(Eneutron[i]>8 && Eneutron[i]<8.5) {
//			cout<<time[i]<<"\t-->\t"<<time_for_neutron<<"\tDist: (m)"<<distance_in_gas<<endl;
//			cout<<Eneutron[i]<<"\t-->\t"<<Eneutronnew[i]<<endl<<endl;
//		}
		if(time[i]<=0) {
			time[i]=0;
			Eneutron[i]=0;
		}
//
		for(int m=0;m<10;m++) {
			Edeuteron=1.0335*Eneutronnew[i]-3.5134;//From zero degrees scattering equation
			fraction=1.-(Edeuteron-5.)/0.21;
			distance_in_gas=fraction*7.98*1e-2;
			neutron_distance=distance_havar-distance_in_gas;
			time_for_beam=distance_in_gas/sqrt((2.*Edeuteron*1.6e-13)/(1.67e-27));//t=d/v
			neutron_cell_time=distance_in_gas/sqrt((2.*Eneutron[i]*1.6e-13)/(1.67e-27));
			time_for_neutron=time[i]-time_for_beam;
			Eneutronnew[i]=(0.5*1.67e-27*TMath::Power(neutron_distance,2.)/(time_for_neutron*time_for_neutron))/1.6e-13;
		}
//
	}
//	TCanvas *c3 = new TCanvas("c3","NE213",600,0,800,1200);
	TCanvas *c3 = new TCanvas();
	c3->cd();
	TGraph *gEn = new TGraph(nbins,Eneutron,neutron);
	gEn->SetTitle(";Neutron energy (MeV);Counts");
	gEn->GetXaxis()->SetRangeUser(6.,14.);
	gEn->SetLineColor(kRed);

//	gEn->Draw();
//	TGraph *gEnnew = new TGraph(nbins,Eneutronnew,neutron);
	cout<<"Here "<<nbins<<endl;
	double errneutron[4096]={0.};
	for(int i=0;i<nbins;i++) errneutron[i]=sqrt(neutron[i]);
	TGraphErrors *gEnnew = new TGraphErrors(nbins,Eneutronnew,neutron,errEneutronnew,errneutron);
//	TGraphErrors *gEnnew = new TGraphErrors(nbins,Eneutron,neutron,errEneutronnew,errneutron);
//	TGraphErrors *gEnnew = new TGraphErrors(nbins,Eneutronnew,neutron,errEneutronnew,errneutron);
	gEnnew->SetTitle(";Neutron energy (MeV);Counts");
	gEnnew->GetXaxis()->SetLimits(3.,15);
	gEnnew->GetXaxis()->SetRangeUser(8.,8.6);
	gEnnew->SetMaximum(250);
	gEnnew->SetLineColor(kMagenta);
	gEnnew->Draw();
	TGraphSmooth *gs = new TGraphSmooth("normal");
	TGraph *grout=gs->SmoothKern(gEn,"normal",0.02);
	grout->SetLineColor(kRed);
	grout->SetLineWidth(3);
//	grout->Draw("SAME");
//
/*
	TGraphSmooth *gsnew = new TGraphSmooth("normal");
	TGraph *groutnew=gsnew->SmoothKern(gEnnew,"normal",0.02);
	groutnew->SetLineColor(kMagenta);
	groutnew->SetLineWidth(3);
//	groutnew->Draw("SAME");
*/
/*	TSpline5 *spline = new TSpline5("spline",gEn->GetX(),gEn->GetY(),gEn->GetN());
	spline->SetLineColor(kRed);
	spline->Draw("lcsame");
*/

//	hEn->GetXaxis()->SetRange(8,10.);
	ifstream in_th;
//	in_th.open("new5.36MeV.txt");
//	in_th.open("dsk303_0.06pc");
//	in_th.open("dsk303");
//	in_th.open("ddn5p362021.txt");
	in_th.open("ddn5p362021_p03air.txt");
	TH1F *hth = new TH1F("hth","",1000,5.,10.);
	TRandom3 *rndm3 = new TRandom3();
	while(in_th.good()) {
		double En,theta;
		in_th>>En>>theta;
		if(theta<0.05) {
//			En=rndm3->Gaus(En,0.004*En);
//			En=rndm3->Gaus(En,0.002*En);
			double vn=sqrt((2.*En*1.6e-13)/(2*1.67e-27));
			for(int l=0;l<10;l++) {
				double depth=(-1+2.*rndm3->Rndm())*2.54*1e-2;
				double timeadjust=(distance_havar-0.0399+depth)/vn;//ns
				double Eadjust=((0.5*2*1.67e-27*(distance_havar-0.0399)*(distance_havar-0.0399))/(timeadjust*timeadjust))/1.6e-13;
//				hth->Fill(Eadjust,1.35*0.0085*2);
				hth->Fill(Eadjust,1.449*0.0085*2);
//				hth->Fill(En,0.085*2);
//				cout<<depth<<"\t"<<En<<"\t"<<Eadjust<<endl;
			}
		}
	}
//	hth->Scale(0.1);
	hth->Draw("SAMEHISTPL");
	hth->SetLineWidth(3);
	hth->GetXaxis()->SetRangeUser(8,8.6);
//

	TLegend *leg2= new TLegend();
	leg2->AddEntry(hth,"GEANT4");
//	leg2->AddEntry(gEn,"Raw NE213");
	leg2->AddEntry(gEnnew,"Corrected NE213");
	leg2->Draw("SAME");
	gPad->SetTicks(1,1);
	return;
}

