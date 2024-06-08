#include <TFile.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TStyle.h>
#include <TLine.h>
#include <TText.h>
#include <TPaletteAxis.h>
#include <iostream>

void modify_canvas(const char* input_filename, const char* output_filename) {
    // Open the input ROOT file
    TFile *inputFile = TFile::Open(input_filename);

    // Check if the input file is opened successfully
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << input_filename << std::endl;
        return;
    }

    // Retrieve the canvas named "Y"
    TCanvas *canvas = (TCanvas*)inputFile->Get("Y");

    // Check if the canvas is retrieved successfully
    if (!canvas) {
        std::cerr << "Error: Cannot find canvas named 'Y' in the input file" << std::endl;
        inputFile->Close();
        return;
    }

    // Draw the canvas
    canvas->cd();

    // Assuming the canvas contains a histogram named "EffMap"
    TH2 *hist = (TH2*)canvas->GetPrimitive("EffMap"); // Replace "EffMap" with the actual name of the histogram if different

    if (!hist) {
        std::cerr << "Error: Cannot find histogram named 'EffMap' in the canvas" << std::endl;
        inputFile->Close();
        return;
    }

    // Remove the statistics box
    hist->SetStats(0);

    // Modify the histogram (title, axis labels, etc.)
    hist->SetTitle("Double tracks efficiency");
    hist->GetXaxis()->SetTitle("ID neuron");
    hist->GetYaxis()->SetTitle("Particle's class");

    // Customize the y-axis labels
    hist->GetYaxis()->SetBinLabel(1, "1 GeV");
    hist->GetYaxis()->SetBinLabel(2, "3 GeV");
    hist->GetYaxis()->SetBinLabel(3, "10 GeV");
    hist->GetYaxis()->SetBinLabel(4, "1 GeV - 1 GeV");
    hist->GetYaxis()->SetBinLabel(5, "1 GeV - 3 GeV");
    hist->GetYaxis()->SetBinLabel(6, "1 GeV - 10 GeV");
    hist->GetYaxis()->SetBinLabel(7, "3 GeV - 3 GeV");
    hist->GetYaxis()->SetBinLabel(8, "3 GeV - 10 GeV");
    hist->GetYaxis()->SetBinLabel(9, "10 GeV - 10 GeV");

    // Draw the histogram with a color palette
    gStyle->SetPalette(kRainBow); // or any other color palette
    gStyle->SetOptStat(0); // Ensure no statistics box is shown
    hist->Draw("COLZ");

    // Set the title for the color bar
    canvas->Update(); // Update the canvas to access the palette
    TPaletteAxis *palette = (TPaletteAxis*)hist->GetListOfFunctions()->FindObject("palette");
    if (palette) {
        palette->SetTitle("Efficiency");
        palette->SetTitleOffset(1.2); // Adjust the offset to ensure visibility
        palette->SetLabelSize(0.03); // Adjust label size if necessary
        palette->GetAxis()->SetTitleSize(0.03); // Adjust title size if necessary
    }

    // Draw a vertical dashed line separating x=9 and x=10 bins
    double x_value = hist->GetXaxis()->GetBinUpEdge(10); // Get the upper edge of the bin 9
    TLine *line = new TLine(x_value, hist->GetYaxis()->GetXmin(), x_value, hist->GetYaxis()->GetXmax());
    line->SetLineStyle(2); // Dashed line
    line->SetLineColor(kBlack); // Black color
    line->SetLineWidth(2); // Set line thickness
    line->Draw();

    // Check for bins with value > 0.8 and add text for all bins in the same column, excluding 0 values
    int nBinsX = hist->GetNbinsX();
    int nBinsY = hist->GetNbinsY();
    for (int i = 1; i <= nBinsX; ++i) {
        bool highlightColumn = false;
        for (int j = 1; j <= nBinsY; ++j) {
            double binContent = hist->GetBinContent(i, j);
            if (binContent > 0.8) {
                highlightColumn = true;
                break;
            }
        }
        if (highlightColumn) {
            for (int j = 1; j <= nBinsY; ++j) {
                double binContent = hist->GetBinContent(i, j);
                if (binContent > 0) {
                    double x = hist->GetXaxis()->GetBinCenter(i);
                    double y = hist->GetYaxis()->GetBinCenter(j);
                    int value = static_cast<int>(binContent * 100);
                    TText *text = new TText(x, y, Form("%d", value));
                    text->SetTextAlign(22); // Center align
                    text->SetTextSize(0.03); // Adjust text size as needed
                    text->Draw();
                }
            }
        }
    }

    // Update the canvas
    canvas->Modified(); // Mark the canvas as modified
    canvas->Update(); // Ensure all changes are visible

    // Open the output ROOT file in update mode
    TFile *outputFile = TFile::Open(output_filename, "UPDATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Unable to open the output ROOT file." << std::endl;
        inputFile->Close();
        return;
    }

    // Write the updated canvas to the output file
    canvas->Write("Y_double_tracks", TObject::kOverwrite);
    outputFile->Close();

    // Close the input file
    inputFile->Close();
}

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input ROOT file> <output ROOT file>" << std::endl;
        return 1;
    }
    modify_canvas(argv[1], argv[2]);
    return 0;
}
//root -l -b -q 'modify_canvas.C("Risultati_singola_traccia.root", "modified_canvas.root")'
