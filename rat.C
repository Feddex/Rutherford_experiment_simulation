#include <iostream>
#include <fstream>
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TROOT.h"

void setStyle(){
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(57);
    gStyle->SetOptTitle(0);
    gStyle->SetOptFit(1); // Enable display of fit information
}

double fitFunction(double* x, double* par) {
    double b = x[0];
    double n = par[0];
    double sigma = par[1];
    if (sin(b / 2) == 0) return 0; // Avoid division by zero
    double pdf = n / (2 * sigma * sin(b / 2) * sin(b / 2));
    return pdf;
}

void Ang() {
    setStyle();

    // Create canvas
    TCanvas *A = new TCanvas("A", "Deflection angular currencies' distribution");

    // Create histogram
    TH1F *h1 = new TH1F("h1", "Angulars", 1000, -0.05, 2 * 3.14159265359 + 0.05);

    // Read data from file and fill histogram
    std::ifstream in("ang10k.txt");
    Float_t y;
    while (1) {
        in >> y;
        if (!in.good()) break;
        h1->Fill(y);
    }
    in.close();

    // Histogram styling
    h1->GetXaxis()->SetTitle("Angle of deflection");
    h1->GetYaxis()->SetTitleOffset(1.);
    h1->GetYaxis()->SetTitle("Currencies");
    h1->SetFillColor(kBlue);
    h1->SetLineWidth(0.1);

    // Define the fit function
    TF1 *fitFunc = new TF1("fitFunc", fitFunction, -0.001, 2 * 3.14159265359+0.001, 2);
    fitFunc->SetParameters(2.25E-25, 1.44E-25); // Adjust initial parameter values for n and sigma
    fitFunc->SetParNames("n", "sigma");
    fitFunc->SetLineColor(kRed); // Set the fit line color to red
    fitFunc->SetLineWidth(2); // Set the fit line width to 1

    // Perform the fit
    h1->Fit("fitFunc", "R");

    // Draw the histogram
    h1->Draw();

    // Customize and display fit results on the plot
    gStyle->SetOptFit(1111); // Display all fit parameters and statistics

    // Update canvas to show fit results
    A->Update();

    // Get the TPaveStats object and move it to a better position
    TPaveStats *stats = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
    if (stats) {
        stats->SetX1NDC(0.7); // X position of the stats box
        stats->SetX2NDC(0.9); // Width of the stats box
        stats->SetY1NDC(0.7); // Y position of the stats box
        stats->SetY2NDC(0.9); // Height of the stats box
    }
    A->Modified();
    A->Update();
}
