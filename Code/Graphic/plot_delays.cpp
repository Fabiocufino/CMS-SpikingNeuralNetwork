#include <TFile.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TLine.h>
#include <iostream>

int plot_delays() {
    // Define the names of the ROOT files
    const char* input_file_name = "Risultati_singola_traccia.root";
    const char* output_file_name = "modified_canvas.root";

    // Open the input ROOT file
    TFile* inputFile = TFile::Open(input_file_name);

    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Unable to open the input ROOT file." << std::endl;
        return 1;
    }

    // Retrieve the TCanvas
    TCanvas* c = (TCanvas*)inputFile->Get("Delta Delay");
    if (!c) {
        std::cerr << "Error: TCanvas 'Delta Delay' not found in the input file." << std::endl;
        return 1;
    }

    // Retrieve the TH2D histogram
    TH2D* hist = (TH2D*)c->GetPrimitive("Delta Delay Histogram");
    if (!hist) {
        std::cerr << "Error: TH2D histogram 'Delta Delay Histogram' not found in the canvas." << std::endl;
        return 1;
    }

    // Multiply all entries by 10^9
    hist->Scale(1e9);

    // Modify the histogram
    hist->GetXaxis()->SetTitle("ID neuron");
    hist->GetYaxis()->SetTitle("ID Tracking Layer");
    hist->SetTitle("Delta Delays [ns]");

    // Set the range for the y-axis to plot only bins y <= 10
    hist->GetYaxis()->SetRangeUser(0, 10);

    // Draw the vertical line at the edge between bins 9 and 10
    double x_line_position = hist->GetXaxis()->GetBinUpEdge(9);
    TLine* line = new TLine(x_line_position, 0, x_line_position, 10);
    line->SetLineStyle(2);  // Dashed line
    line->SetLineWidth(3);  // Thick line
    line->Draw();

    gStyle->SetPalette(kRainBow); // or any other color palette
    gStyle->SetOptStat(0); // Ensure no statistics box is shown

    // Update and draw the canvas
    c->Update();
    c->Draw();

    // Open the output ROOT file in update mode
    TFile* outputFile = TFile::Open(output_file_name, "UPDATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Unable to open the output ROOT file." << std::endl;
        return 1;
    }

    // Write the updated canvas to the output file
    c->Write("Updated_Delta_Delay", TObject::kOverwrite);
    outputFile->Close();

    // Close the input file
    inputFile->Close();

    return 0;
}
