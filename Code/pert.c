#include <boost/random.hpp>
#include <boost/random/beta_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include "TH1F.h"
#include "TCanvas.h"

using namespace boost::random;
using namespace std;

mt19937 rng;

int main(int argc, char const *argv[])
{
    float a = 0;
    float b = 1;
    float c = 1;

    float alpha = (4 * b + c - 5 * a) / (c - a);
    float beta = (5 * c - a - 4 * b) / (c - a);
    beta_distribution<> betaDistribution(alpha, beta);
    variate_generator<mt19937&, beta_distribution<> > betaGenerator(rng, betaDistribution);

    // Generate random numbers
    const int num_samples = 1000;
    vector<float> random_numbers;
    for (int i = 0; i < num_samples; ++i) {
        float random = betaGenerator();
        random_numbers.push_back(random);
    }

    // Plot the histogram using ROOT
    TH1F* hist = new TH1F("hist", "Beta Distribution", 100, 0, 1);
    for (const auto& num : random_numbers) {
        hist->Fill(num);
    }

    // Create a canvas and draw the histogram
    TCanvas* canvas = new TCanvas("canvas", "Canvas", 800, 600);
    hist->Draw();

    // Save the canvas to a PDF file
    canvas->SaveAs("output.pdf");

    return 0;
}
