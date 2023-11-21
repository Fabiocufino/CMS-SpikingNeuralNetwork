#include <iostream>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TApplication.h>
#include <TBrowser.h>

int main(int argc, char **argv) {
    TApplication theApp("App", &argc, argv);

    // Create some example data
    const int numPoints = 100;
    double x[numPoints];
    double y1[numPoints], y2[numPoints], y3[numPoints], y4[numPoints], y5[numPoints], y6[numPoints];

    for (int i = 0; i < numPoints; ++i) {
        x[i] = i;
        y1[i] = sin(0.1 * x[i]);
        y2[i] = cos(0.1 * x[i]);
        y3[i] = 0.5 * sin(0.2 * x[i]);
        y4[i] = 0.5 * cos(0.2 * x[i]);
        y5[i] = 0.2 * sin(0.3 * x[i]);
        y6[i] = 0.2 * cos(0.3 * x[i]);
    }
    TCanvas *canvas = new TCanvas("canvas", "Multiple Graphs", 800, 600);

    // Create TGraph objects for each function
    TGraph *graph1 = new TGraph(numPoints, x, y1);
    TGraph *graph2 = new TGraph(numPoints, x, y2);
    TGraph *graph3 = new TGraph(numPoints, x, y3);
    TGraph *graph4 = new TGraph(numPoints, x, y4);
    TGraph *graph5 = new TGraph(numPoints, x, y5);
    TGraph *graph6 = new TGraph(numPoints, x, y6);

    // Create a TMultiGraph to hold multiple graphs
    TMultiGraph *mg = new TMultiGraph();
    mg->Add(graph1, "l");
    mg->Add(graph2, "l");
    mg->Add(graph3, "l");
    mg->Add(graph4, "l");
    mg->Add(graph5, "l");
    mg->Add(graph6, "l");

    // Create a TCanvas to display the graphs

    // Create a TBrowser with check boxes for each graph
    // TBrowser *browser = new TBrowser("Browser");

    // Set the TMultiGraph as the object to be displayed
    // browser->SetObject(mg);
    mg -> Draw("a");

    // Event loop to handle user input
    theApp.Run();

    return 0;
}
