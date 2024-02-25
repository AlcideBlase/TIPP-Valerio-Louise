/*
Protons File : File HEPData-1569102768-v1.root
Table 5 : Pb-Pb collision 
Table 6 : pp collision

Mesons File : File HEPData-ins1762368-v1.root
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
#include "TLatex.h"

//-----------------------------------------------------------------------------------------------------------------------------------------------
//Initialisation of differents variables
TString file = Form("./data/HEPData-ins1762368-v1.root"); //You can choose the file here
TString table = Form("Table 4"); //You can choose the table here
const int n_histogram = 1; //You can choose the number of histograms

//Variables that will depends on the file and the table taked
TString title = Form(" ");
TString expo_fit_para = Form(" ");
TString boltz_fit_para = Form(" ");
TString power_fit_para = Form(" ");
TString levy_fit_para = Form(" ");
TString save_name = Form("./output/ParticulesFit");


//Definition of the differents masses in Gev
double masse_phi = 1.0194; 
double masse_proton = 0.93827208816;
double masse = masse_proton;


//List for the value of the integral, a list is needed if there is more than 1 histogramm.
double integral_expo[n_histogram];
double integral_boltz[n_histogram];
double integral_law[n_histogram];
double integral_levy[n_histogram];


// We add the errors
double err1Value = 0;
double err2Value = 0;
double combinedError = 0;
double binError = 0;


double xmin = 0;
double xmax = 0;


//Parameter for the legends
double dNdy_levy = 0;
double temp_levy = 0;
double n_levy = 0;
double error_dN = 0;
double error_T = 0;
double error_n = 0;


//Definition of the differents function
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
	double A = par[0]; 
	double p0 = par[1];
	double n = par[2];
	double B = 1+pow((x[0]/p0),2);
	return A*x[0]/(pow(B,n));
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


//Suppresion of the bin without value that are egal to 0
//We think that if there is less value in one of the table it will appear as a value exactly egal to 0 at the end.
int ghost_data(TH1F *Yields){
		double NoNull = 0;
		double step = 0;
		int nmax = 0;
		while(NoNull == 0) {
			NoNull += abs(Yields->GetBinContent(Yields->GetNbinsX()-step));
			step++;
		}
	step--;
	nmax = Yields->GetNbinsX()-step;
	return nmax;
}

double integrale_histo(TH1F *Yield, int number_bins){
	double integrale = 0;
	double largeur = 0;
	double aire= 0;
	double hauteur=0;
	TAxis *xAxis = Yield->GetXaxis();
	for (int k = 1; k <= number_bins; ++k) {
		largeur = abs(xAxis->GetBinCenter(k)-xAxis->GetBinCenter(k-1));
		hauteur = Yield->GetBinContent(k);
		aire = largeur * hauteur;
		integrale += aire;
	}
	return integrale;
}

void ParticulesFit(){
	//We recover the file and the table defined before
	TFile *myFile = new TFile(file);
	TDirectoryFile *MesonDirectoryFile = (TDirectoryFile*)myFile->Get(table);


	if (table == "Table 3" && file == "./data/HEPData-ins1762368-v1.root") { 
		expo_fit_para = "IE+"; //The best fit method depends on the function and on the file so we need to take that into account
		boltz_fit_para = "IE+";
		power_fit_para = "IE+";
		levy_fit_para = "IE+";
		title = "p_{T} distributions of #phi meson measured in Pb-Pb collisions at #sqrt{s}= 5.02 TeV";
		masse = masse_phi; //The mass need to be different if we study proton of phi mesons
		save_name = save_name+"_meson_table_3.pdf";
	}
		
	else if(table == "Table 4" && file == "./data/HEPData-ins1762368-v1.root") {
		expo_fit_para = "IEM+";
		boltz_fit_para = "IEM+";
		power_fit_para = "EM+";
		levy_fit_para = "IEM+";
		title = "p_{T} distributions of #phi meson measured in pp collisions at #sqrt{s}= 5.02 TeV";
		masse = masse_phi;
		save_name = save_name+"_meson_table_4.pdf";
	}
		
	else if(table == "Table 5" && file == "./data/HEPData-1569102768-v1.root") {
		expo_fit_para = "REM+";
		boltz_fit_para = "REM+";
		power_fit_para = "REM+";
		levy_fit_para = "IE+";
		title = "p_{T} distributions of p-#bar{p} measured in Pb-Pb collisions at #sqrt{s}= 5.02 TeV";
		masse = masse_proton;
		save_name = save_name+"_proton_table_5.pdf";
	}
	
	else if(table == "Table 6" && file == "./data/HEPData-1569102768-v1.root") {
		expo_fit_para = "IEM+";
		boltz_fit_para = "IEM+";
		power_fit_para = "EM+";
		levy_fit_para = "EM+";
		title = "p_{T} distributions of p-#bar{p} measured in pp collisions at #sqrt{s}= 5.02 TeV";
		masse = masse_proton;
		save_name = save_name+"_proton_table_6.pdf";
	}

	else {
		expo_fit_para = "IEM+";
		boltz_fit_para = "IEM+";
		power_fit_para = "IEM+";
		levy_fit_para = "IEM+";
		title = "p_{T} distributions of a collision";
		save_name = save_name+table;
		save_name = save_name+".pdf";
	}


	//Creation of the Canvas
	TCanvas *firstCanvas = new TCanvas("firstCanvas","Histogrammes ",1200,800);
	firstCanvas->SetMargin(0.15,0.05,0.15,0.05);
	firstCanvas->cd();
	firstCanvas->SetLogy();
	TH2F *blankHisto = new TH2F("histogram",title,21,0,20,1000,1e-7,100);
	blankHisto->SetXTitle("p_{T} (GeV/c)");
	blankHisto->GetXaxis()->SetTitleSize(0.04);
	blankHisto->GetXaxis()->SetTitleOffset(1.2);
	blankHisto->SetYTitle("#frac{d^{2}N}{dydp_{T}} (GeV/c)^{-1} ");
	blankHisto->GetYaxis()->SetTitleSize(0.04);
	blankHisto->GetYaxis()->SetTitleOffset(1.6);
	blankHisto->SetTitleOffset(1.5);
	blankHisto->SetStats(0);
	blankHisto->Draw();


	for(int i = 1; i <= n_histogram; ++i) { //We loop for all the histograms in the table

		//Name of the histogram and their errors
		TString histName = Form("Hist1D_y%d", i);
		TString histErrName1 = Form("Hist1D_y%d_e1", i);
		TString histErrName2 = Form("Hist1D_y%d_e2", i);

		// We get the histogram and their errors
		TH1F *myYields = (TH1F*)MesonDirectoryFile->Get(histName.Data());
		TH1F *myErr1 = (TH1F*)MesonDirectoryFile->Get(histErrName1.Data());
		TH1F *myErr2 = (TH1F*)MesonDirectoryFile->Get(histErrName2.Data());
		
		blankHisto->GetYaxis()->SetRangeUser(1e-7,myYields->GetBinContent(1)*10); //We define here the range of the Y axix that will depends on the value of the histogram.

		if (table == "Table 4" && file == "./data/HEPData-ins1762368-v1.root") { // The table 4 of the Meson file got some corrupted data so we change them manually we use the data on the online version of the file
			myYields->SetBinContent(5,0.0167252);
			myErr1->SetBinContent(5,0.000210925);
			myErr2->SetBinContent(5,0.000967905);
			myYields->SetBinContent(6,0.015271);
			myErr1->SetBinContent(6,0.0001848);
			myErr2->SetBinContent(6,0.000803366);
		}

		int nmax = ghost_data(myYields); //We erase the data at the end that are not really there
		myYields->GetXaxis()->SetRange(0,nmax);


		//Set the different range for the fit and the integrals
		TAxis *xaxis = myYields->GetXaxis(); //Take the Xaxis
		xmin = xaxis->GetBinCenter(1); //Choose the minimum for the x axis by choosing the x value of the first bin (the 0 bin isn't a real value)
		xmax = xaxis->GetBinCenter(nmax); //Choose the maximum value fot the x axis by using the while loop who verify if there is dead value at the end.


		// Add the errors
		for (int iBin = 1; iBin <= myYields->GetNbinsX(); ++iBin) { //We start at 1 because the first bin isn't a real one.
			err1Value = myErr1->GetBinContent(iBin);
			err2Value = myErr2->GetBinContent(iBin);
			combinedError = sqrt(err1Value * err1Value + err2Value * err2Value); //That's how we add the statistical and systematic error : by using quadrature sommation
			myYields->SetBinError(iBin, combinedError);
		}
		
		myYields->SetLineColor(i); // Change the color of every histogramm
		
		myYields->SetStats(0);
		myYields->Draw("SAME"); // "SAME" to write on the same canvas everything.

		// Create function for every histogram
		TF1 * func_levy = new TF1(Form("func_levy_%d", i),levy,0,20,4);
		TF1 * func_expo = new TF1(Form("func_expo_%d", i),exponentielle,0,20,3);
		TF1 * func_boltz = new TF1(Form("func_boltz_%d", i),boltzmann,0,20,3);
		TF1 * func_law = new TF1(Form("func_law_%d", i),power_law,0,20,3);


		// Write the parameters of the function for the fit
		//The parameters are define by a first run done on other data and with a step really low so that the parameter are set on value that are on a range physically appropriate
		func_expo->SetParameter(0,1.0);
		func_expo->SetParameter(1,1.0);
		func_expo->SetParameter(2,masse);
		func_expo->SetParLimits(2,masse,masse); //We fixed the value of the mass by setting the below and above limit as it value

		func_boltz->SetParameter(0,1.0);
		func_boltz->SetParameter(1,1.0);
		func_boltz->SetParameter(2,masse);
		func_boltz->SetParLimits(2,masse,masse);

		func_law->SetParameter(0,1.0);
		func_law->SetParameter(1,1.0);
		func_law->SetParameter(2,1.0);

		func_levy->SetParameter(0,0.5);
		func_levy->SetParameter(1,7.0);
		func_levy->SetParameter(2,2.0);
		func_levy->SetParameter(3,masse);
		func_levy->SetParLimits(3,masse,masse);

		func_expo->SetLineColor(1); //We change the color of every fit
		func_boltz->SetLineColor(2);
		func_law->SetLineColor(3);
		func_levy->SetLineColor(4);


		//ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2"); //We set the algorithm for the minimalization


		//We fit the function
		//Fitting two times made sometimes better results
		myYields->Fit(func_expo, expo_fit_para+"Q", "", xmin, xmax);  //We add "Q" for the first run so that he will be quiet and wont appear in the end
		myYields->Fit(func_expo, expo_fit_para, "", xmin, xmax); 
		
		myYields->Fit(func_boltz, boltz_fit_para+"Q", "", xmin, xmax);
		myYields->Fit(func_boltz, boltz_fit_para, "", xmin, xmax);
		
		myYields->Fit(func_law, power_fit_para+"Q", "", xmin, xmax);
		myYields->Fit(func_law, power_fit_para, "", xmin, xmax);

		myYields->Fit(func_levy, levy_fit_para+"Q", "", xmin, xmax);
		myYields->Fit(func_levy, levy_fit_para, "", xmin, xmax);


		//Calcul of the differents integrals
		integral_expo[i] = func_expo->Integral(0,20);
		integral_boltz[i] = func_boltz->Integral(0,20);
		integral_law[i] = func_law->Integral(0,20);
		integral_levy[i] = func_levy->Integral(0,20);
		
		
		//We draw the function on all the interval 
		func_expo->Draw("same");
		func_boltz->Draw("same");
		func_law->Draw("same");
		func_levy->Draw("same");


		//We display the value of the differents integrals
		std::cout << "Integral of the expo function for the histogramm number : " << i << " : " << integral_expo[i] << std::endl;
		std::cout << "Integral of the boltzmann function for the histogramm number : " << i << " : " << integral_boltz[i] << std::endl;
		std::cout << "Integral of the power_law function for the histogramm number : " << i << " : " << integral_law[i] << std::endl;	
		std::cout << "Integral of the Levy function for the histogramm number : " << i << " : " << integral_levy[i] << std::endl;
		std::cout << "Integrale of the histogram : " << integrale_histo(myYields, nmax) << std::endl;
		
		// Legends
		TLegend *legend = new TLegend(0.6,0.6,0.95,0.95);

		//Parameters and their error (of the Levy-Tsallis function for example)
		double dNdy_levy = func_levy->GetParameter(0); // dN/dy parameter 
		double temp_levy = func_levy->GetParameter(2); // temperature parameter
		double n_levy = func_levy->GetParameter(1); // n parameter
		double error_dN = func_levy->GetParError(0); // error on dN/dy
		double error_T = func_levy->GetParError(2); //error on temperature
		double error_n = func_levy->GetParError(1); //error on n
		double temp_levy_multiplied = temp_levy * 1000; //the temperature is multiplied by X1000
		double error_T_multiplied = error_T * 1000; //error of the temperature X1000

		//Print the parameter +/- their error
		Char_t chardNdy[45]; //Choose the size of the parameter
		snprintf(chardNdy,45,"dN/dy_{Levy} = %.5f #pm %.5f", dNdy_levy, error_dN);
		Char_t charTemp[45];
		snprintf(charTemp, 45, "T_{Levy}   = %.0f #pm %.0f MeV", temp_levy_multiplied, error_T_multiplied); 
		Char_t charn[45]; 
		snprintf(charn,45,"n_{Levy}  = %.1f #pm %.1f ",n_levy, error_n);

		//Draw the parameters in LaTeX
		TLatex *fitText = new TLatex(0.75,0.52,chardNdy);
		fitText->SetNDC(1);
		fitText->SetTextSize(0.03);
		fitText->DrawLatex(0.65,0.65,chardNdy);
		fitText->DrawLatex(0.65,0.60,charTemp);
		fitText->DrawLatex(0.65,0.55,charn);
  	
		TLegend *tleg = new TLegend(0.65,0.7,0.90,0.8);
		tleg->AddEntry(func_expo,"Exponential fit","l");
		tleg->AddEntry((TObject*)nullptr, "", ""); //We put a empty legend to make space between legends
		tleg->AddEntry(func_boltz,"Boltzmann fit","l");
		tleg->AddEntry((TObject*)nullptr, "", "");
		tleg->AddEntry(func_law,"Power law fit","l");
		tleg->AddEntry((TObject*)nullptr, "", "");
		tleg->AddEntry(func_levy,"Levy fit","l");
		
		tleg->SetTextSize(0.03);
		tleg->SetBorderSize(0);
		tleg->Draw();
	}
	firstCanvas->Draw(); //We display the Canvas
	firstCanvas->SaveAs(save_name); //We save the canva as a pdf image
}
