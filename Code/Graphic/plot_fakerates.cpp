#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TLine.h>
#include <iostream>
#include <algorithm>

int plot_fakerates() {
    // Define the number of histograms
    const int n_histograms = 20;

    // Define the names of the ROOT files
    const char* file1_name = "Risultati_singola_traccia.root";
    const char* file2_name = "Risultati_doppia_traccia.root";
    const char* output_file_name = "modified_canvas.root";

    // Open the ROOT files
    TFile* file1 = TFile::Open(file1_name);
    TFile* file2 = TFile::Open(file2_name);

    if (!file1 || !file2) {
        std::cerr << "Error: Unable to open one or both ROOT files." << std::endl;
        return 1;
    }

    // Arrays to store the fake rates
    double fake_rates_file1[n_histograms];
    double fake_rates_file2[n_histograms];

    // Retrieve the histograms and get the fake rates
    for (int i = 0; i < n_histograms; ++i) {
        std::string hist_name = "FakeRate" + std::to_string(i);

        TH1F* hist_file1 = (TH1F*)file1->Get(hist_name.c_str());
        TH1F* hist_file2 = (TH1F*)file2->Get(hist_name.c_str());

        if (!hist_file1 || !hist_file2) {
            std::cerr << "Error: Histogram " << hist_name << " not found in one of the files." << std::endl;
            return 1;
        }

        fake_rates_file1[i] = hist_file1->GetBinContent(1);
        fake_rates_file2[i] = hist_file2->GetBinContent(1);
    }

    // Close the ROOT files
    file1->Close();
    file2->Close();

    // Create a plot
    TCanvas* c = new TCanvas("FakeRates", "Fake Rates", 800, 600);
    TMultiGraph* mg = new TMultiGraph();

    // Convert the fake rates to ROOT TGraph
    TGraph* graph_file1 = new TGraph(n_histograms);
    TGraph* graph_file2 = new TGraph(n_histograms);

    for (int i = 0; i < n_histograms; ++i) {
        graph_file1->SetPoint(i, i, fake_rates_file1[i]);
        graph_file2->SetPoint(i, i, fake_rates_file2[i]);
    }

    // Set marker styles and colors
    graph_file1->SetMarkerStyle(20);
    graph_file1->SetMarkerColor(kBlue);
    graph_file1->SetLineColor(kBlue);

    graph_file2->SetMarkerStyle(21);
    graph_file2->SetMarkerColor(kRed);
    graph_file2->SetLineColor(kRed);

    // Add graphs to the multigraph
    mg->Add(graph_file1);
    mg->Add(graph_file2);

    // Draw the multigraph
    mg->Draw("APL");
    mg->SetTitle("Fake Rates;Neuron;Fake Rate");
    mg->GetXaxis()->SetLimits(-1, n_histograms);
    mg->GetYaxis()->SetRangeUser(0, std::max(*std::max_element(fake_rates_file1, fake_rates_file1 + n_histograms), *std::max_element(fake_rates_file2, fake_rates_file2 + n_histograms)) * 1.1);

    // Create a legend
    TLegend* legend = new TLegend(0.1, 0.7, 0.3, 0.9);
    legend->AddEntry(graph_file1, "Single Tracks", "lp");
    legend->AddEntry(graph_file2, "Double Tracks", "lp");
    legend->Draw();

    // Draw the vertical dashed thick line at x = 9.5
    TLine* line = new TLine(9.5, 0, 9.5, mg->GetYaxis()->GetXmax());
    line->SetLineStyle(2);  // Dashed line
    line->SetLineWidth(3);  // Thick line
    line->Draw();

    // Draw the canvas
    c->Update();
    c->Draw();

    // Open the existing ROOT file in update mode
    TFile* outputFile = new TFile(output_file_name, "UPDATE");
    if (!outputFile) {
        std::cerr << "Error: Unable to open the output ROOT file." << std::endl;
        return 1;
    }

    // Write the updated canvas to the file
    c->Write("", TObject::kOverwrite);
    outputFile->Close();

    return 0;
}
