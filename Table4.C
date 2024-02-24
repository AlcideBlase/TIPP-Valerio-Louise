/*
Protons : File HEPData-1569102768-v1.root
Table 5 : Pb-Pb collision 
Table 6 : pp collision

Mesons: File HEPData-ins1762368-v1.root
Table 3 : Pb-Pb collision
Table 4 : pp  collision
*/

#include <iostream>
#include <cmath>
#include "TString.h"
#include "Math/MinimizerOptions.h"
#include "TLegend.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TDirectoryFile.h"
#include "TCanvas.h"
#include "TF1.h"

TString file = Form("HEPData-ins1762368-v1.root"); //You can choose the file here
TString table = Form("Table 4"); //You can choose the table here

TString titre = Form("");
TString expo_fit_para = Form("");
TString boltz_fit_para = Form("");
TString power_fit_para = Form("");
TString levy_fit_para = Form("");

double masse_phi = 1.0194;

double exponentielle(double *x, double *par) {
	double m=par[2];
	double s1 = sqrt(pow(m,2) +pow(x[0],2));
	double s2 = sqrt(2*pow(m,2) +pow(x[0],2));
	return par[0] * s1 * exp(-s2/ par[1]);
}

double boltzmann(double *x, double *par){
	double m=par[2];
	double s1 = sqrt(pow(m,2) +pow(x[0],2));
	double s2 = sqrt(2*pow(m,2) +pow(x[0],2));
	return par[0] * s1 *s2* exp(-s2/ par[1]);
}

double power_law(double *x, double *par){
	return par[0]*x[0]/(pow((1+pow((x[0]/par[1]),2)),par[2]));
}

double levy(double *x, double *par){
	double masse=par[3];
	double m_T=sqrt(masse*masse +x[0]*x[0]);
	double DN_Dy = par[0];
	double n = par[1];
	double T = par[2];
	double nT= n*T;
	return DN_Dy*((n-1)*(n-2)/(nT*(nT+(n-2)*masse)))*x[0]*pow(1+(m_T-masse)/(nT),-n);
}

void Table4(){
	//We recover the file and the good table (number 3)
	TFile *myFile = new TFile(file);
	TDirectoryFile *MesonDirectoryFile = (TDirectoryFile*)myFile->Get(table);


	if (table == "Table 3" && file == "HEPData-ins1762368-v1.root") { //The best fit method depends on the function and on the file so we need to take that into account
			expo_fit_para = "IE+";
			boltz_fit_para = "IE+";
			power_fit_para = "IE+";
			levy_fit_para = "IE+";
			titre = "p_{T} distributions of phi meson measured in Pb-Pb collisions at sqrt(s)= 5.02 TeV";
		}
		
		else if(table == "Table 4" && file == "HEPData-ins1762368-v1.root") {
			expo_fit_para = "IEM+";
			boltz_fit_para = "IEM+";
			power_fit_para = "EM+";
			levy_fit_para = "IEM+";
			titre = "p_{T} distributions of phi meson measured in pp collisions at sqrt(s)= 5.02 TeV";
		}
		
		else if(table == "Table 5" && file == "HEPData-1569102768-v1.root") {
			expo_fit_para = "REM+";
			boltz_fit_para = "REM+";
			power_fit_para = "REM+";
			levy_fit_para = "IE+";
			titre = "p_{T} distributions of p-pbar measured in Pb-Pb collisions at sqrt(s)= 5.02 TeV";
		}
		
		else if(table == "Table 6" && file == "HEPData-1569102768-v1.root") {
			expo_fit_para = "IEM+";
			boltz_fit_para = "IEM+";
			power_fit_para = "EM+";
			levy_fit_para = "EM+";
			titre = "p_{T} distributions of p-pbar measured in pp collisions at sqrt(s)= 5.02 TeV";
		}
		else {
			expo_fit_para = "IEM+";
			boltz_fit_para = "IEM+";
			power_fit_para = "IEM+";
			levy_fit_para = "IEM+";
			titre = "p_{T} distributions of a collision";
		}

	//Creation of the Canvas
	TCanvas *firstCanvas = new TCanvas("firstCanvas","Histogrammes ",1200,800);
	firstCanvas->SetMargin(0.15,0.05,0.15,0.05);
	firstCanvas->cd();
	firstCanvas->SetLogy();
	TH2F *blankHisto = new TH2F("histogram",titre,21,0,20,1000,1e-7,100);
	blankHisto->SetXTitle("p_{T} (GeV/c)");
	blankHisto->SetStats(0);
	blankHisto->Draw();

	//Creation of the differents variables

	double integral_expo[20]; //List for the value of the integral.
	double integral_boltz[20];
	double integral_law[20];
	double integral_levy[20]; 

	double NoNull = 0; //Variables needed for the suppresion of the bin without value.
	double step = 0;

	// For loop to write et draw the main histogramms et their errors.
	int i = 1;
	//for(int i = 1; i <= 1; ++i) { //You can use the for loop if you want to see on all the histogram
		// Name of the histogram and their errors
		TString histName = Form("Hist1D_y%d", i);
		TString histErrName1 = Form("Hist1D_y%d_e1", i);
		TString histErrName2 = Form("Hist1D_y%d_e2", i);

		// We get the histogram and their errors
		TH1F *myYields = (TH1F*)MesonDirectoryFile->Get(histName.Data());
		TH1F *myErr1 = (TH1F*)MesonDirectoryFile->Get(histErrName1.Data());
		TH1F *myErr2 = (TH1F*)MesonDirectoryFile->Get(histErrName2.Data());
		
		blankHisto->GetYaxis()->SetRangeUser(1e-7,myYields->GetBinContent(1)*10);
		
		if (table == "Table 4" && file == "HEPData-ins1762368-v1.root") { // The table 4 of the Meson file got some corrupted data so we change them manually
			myYields->SetBinContent(5,0.0167252);
			myErr1->SetBinContent(5,0.000210925);
			myErr2->SetBinContent(5,0.000967905);
			myYields->SetBinContent(6,0.015271);
			myErr1->SetBinContent(6,0.0001848);
			myErr2->SetBinContent(6,0.000803366);
		}

		//Suppresion of the bin without value that are egal to 0
		//We think that if there is less value in one of the table it will appear as a value exactly egal to 0 at the end.
		NoNull = 0;
		step = 0;
		while(NoNull == 0) {
			NoNull = NoNull + abs(myYields->GetBinContent(myYields->GetNbinsX()-step));
			step++; 
		}
		step--;  
		int nmax = myYields->GetNbinsX()-step;    
		myYields->GetXaxis()->SetRange(0,nmax);


		//Set the different range for the fit and the integrals
		TAxis *xaxis = myYields->GetXaxis(); //Take the Xaxis
		double xmin = xaxis->GetBinCenter(1); //Choose the minimum for the x axis by choosing the x value of the first bin (the 0 bin isn't a real value)
		double xmax = xaxis->GetBinCenter(nmax); //Choose the maximum value fot the x axis by using the while loop who verify if there is dead value at the end.


		// We add the errors
		double err1Value = 0;
		double err2Value = 0;
		double combinedError = 0;
		double binError = 0;
		for (int iBin = 1; iBin <= myYields->GetNbinsX(); ++iBin) {
			err1Value = myErr1->GetBinContent(iBin); //We put +1 because the first bin isn't a real one.
			err2Value = myErr2->GetBinContent(iBin);
			combinedError = sqrt(err1Value * err1Value + err2Value * err2Value); //That's how we add the statistical and systematic error : by using quadrature
			myYields->SetBinError(iBin, combinedError);
		}

		// Add the errors
		myYields->SetLineColor(i); // Change the color of every histogramm
		
		myYields->SetStats(0);
		myYields->Draw("SAME"); // "same" to write on the same canvas everything.

		// Integral of histogram
		
		double integrale = 0; 
		double largeur = 0;
		double aire= 0;
		double hauteur=0;
		
		for (int k = 1; k <= nmax; ++k) {
			largeur = abs(xaxis->GetBinCenter(k)-xaxis->GetBinCenter(k-1));
    			hauteur = myYields->GetBinContent(k);
    			aire = largeur * hauteur;
    			integrale += aire;
 		}
 		std::cout << "integrale de l'histo : " << integrale << std::endl;
		
		// Create function for every histogram
		TF1 * func_levy = new TF1(Form("func_levy_%d", i),levy,0,20,4);
		TF1 * func_expo = new TF1(Form("func_expo_%d", i),exponentielle,0,20,3);
		TF1 * func_boltz = new TF1(Form("func_boltz_%d", i),boltzmann,0,20,3);
		TF1 * func_law = new TF1(Form("func_law_%d", i),power_law,0,20,3);

		// Write the parameters of the function for the fit
						
		func_expo->SetParameter(0,1.0);
		func_expo->SetParameter(1,1.0);
		func_expo->SetParameter(2,masse_phi);
		func_expo->SetParLimits(2,masse_phi,masse_phi); //We fixed the value of the mass
		
		func_boltz->SetParameter(0,1.0);
		func_boltz->SetParameter(1,1.0);
		func_boltz->SetParameter(2,masse_phi);
		func_boltz->SetParLimits(2,masse_phi,masse_phi);
		
		func_law->SetParameter(0,1.0);
		func_law->SetParameter(1,1.0);
		func_law->SetParameter(2,1.0);
		
		func_levy->SetParameter(0,0.5);
		func_levy->SetParameter(1,7.0);
		func_levy->SetParameter(2,2.0);
		func_levy->SetParameter(3,masse_phi);
		func_levy->SetParLimits(3,masse_phi,masse_phi); 

		//ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2"); //We set the algorithm for the minimalization
		
		func_expo->SetLineColor(1); //We set same color the fit and the histogramm
		func_boltz->SetLineColor(2);
		func_law->SetLineColor(3);
		func_levy->SetLineColor(4);
		
		//We fit the function, RI is for taking the good range with R, and with I to fit with the value of the integral and not the value of the bin
		//Fitting two times made sometimes better results
		myYields->Fit(func_expo, expo_fit_para+"Q", "", xmin, xmax); 
		myYields->Fit(func_expo, expo_fit_para, "", xmin, xmax); 
		
		myYields->Fit(func_boltz, boltz_fit_para+"Q", "", xmin, xmax);
		myYields->Fit(func_boltz, boltz_fit_para, "", xmin, xmax);
		
		myYields->Fit(func_law, power_fit_para+"Q", "", xmin, xmax);
		myYields->Fit(func_law, power_fit_para, "", xmin, xmax);

		myYields->Fit(func_levy, levy_fit_para+"Q", "", xmin, xmax);
		myYields->Fit(func_levy, levy_fit_para, "", xmin, xmax);


		integral_expo[i] = func_expo->Integral(0,20); //Calcul of the integrale
		integral_boltz[i] = func_boltz->Integral(0,20);
		integral_law[i] = func_law->Integral(0,20);
		integral_levy[i] = func_levy->Integral(0,20);
		
		std::cout << "Integral of the expo function for the histogramm number : " << i << " : " << integral_expo[i] << std::endl;//We display the value of the integral
		std::cout << "Integral of the boltzmann function for the histogramm number : " << i << " : " << integral_boltz[i] << std::endl;
		std::cout << "Integral of the power_law function for the histogramm number : " << i << " : " << integral_law[i] << std::endl;	
		std::cout << "Integral of the levy function for the histogramm number : " << i << " : " << integral_levy[i] << std::endl; 
		
		// Légendes
		TLegend *legend = new TLegend(0.6,0.6,0.95,0.95);
		
		double dNdy_levy = func_levy->GetParameter(0); // dN/dy
		double temp_levy = func_levy->GetParameter(2); // temperature
		double n_levy = func_levy->GetParameter(1); // n
		double error_dN = func_levy->GetParError(0);
		double error_T = func_levy->GetParError(2);
		double error_n = func_levy->GetParError(1);
		
		Char_t chardNdy[45];
  		snprintf(chardNdy,45,"dN/dy_{Levy} = %.5f +/- %.5f", dNdy_levy, error_dN);
  		Char_t charTemp[45]; ;
  		snprintf(charTemp,45,"T_{Levy}   = %.2f +/- %.2f (MeV)",temp_levy, error_T);
  		Char_t charn[45]; 
  		snprintf(charn,45,"n_{Levy}  = %.2f +/- %.2f ",n_levy, error_n);

  		TLatex *fitText = new TLatex(0.75,0.52,chardNdy);
  		fitText->SetNDC(1);
  		fitText->SetTextSize(0.03);
  		fitText->DrawLatex(0.65,0.75,chardNdy);
  		fitText->DrawLatex(0.65,0.70,charTemp);
  		fitText->DrawLatex(0.65,0.65,charn);
  	
 		TLegend *tleg = new TLegend(0.65,0.8,0.90,0.9);
 		tleg->AddEntry(func_expo,"Exponential fit","l");
 		tleg->AddEntry(func_boltz,"Boltzmann fit","l");
 		tleg->AddEntry(func_law,"Power law fit","l");
 		tleg->AddEntry(func_levy,"Levy fit","l");
 		tleg->SetTextSize(0.03);
 		//tleg->SetFillColor(10);
  		tleg->SetBorderSize(0);
  		tleg->Draw();
	//}
	firstCanvas->Draw(); //We display the Canvas

}
