#include <TFile.h>
#include <TDirectoryFile.h>
#include <TGraph.h>
#include <TCanvas.h>

void read_HEPData() {
    // Open the ROOT file
    TFile file("HEPData-ins2674351-v1-root.root");
    //Download root file from https://www.hepdata.net/record/ins2674351

    // Check if the file is open
    if (!file.IsOpen()) {
        std::cerr << "Error: Unable to open ROOT file." << std::endl;
        return;
    }

    // Navigate to the TDirectoryFile fig_04
    TDirectoryFile *dir = dynamic_cast<TDirectoryFile*>(file.Get("fig_04"));
    if (!dir) {
        std::cerr << "Error: TDirectoryFile 'fig_04' not found." << std::endl;
        file.Close();
        return;
    }

    // Create a canvas
    TCanvas *canvas = new TCanvas("canvas", "Plot from HEPData", 800, 600);
    canvas->SetLogy(); // Set logarithmic scale for y-axis

    // Loop through each graph
    for (int i = 1; i <= 9; ++i) {
        // Access the TGraph from the TDirectoryFile
        TString graphName = TString::Format("Graph1D_y%d", i);
        TGraph *graph = dynamic_cast<TGraph*>(dir->Get(graphName));
        if (!graph) {
            std::cerr << "Error: TGraph '" << graphName << "' not found in 'fig_04'." << std::endl;
            continue; // Skip to the next graph
        }

        // Draw the graph
        graph->Draw("ALP");

        // Customize the graph appearance
        TString yAxisLabel;
        switch (i) {
            case 1: yAxisLabel = "Mjj"; break;
            case 2: yAxisLabel = "Mjb"; break;
            case 3: yAxisLabel = "Mbb"; break;
            case 4: yAxisLabel = "Mje"; break;
            case 5: yAxisLabel = "Mbe"; break;
            case 6: yAxisLabel = "Mby"; break;
            case 7: yAxisLabel = "Mjm"; break;
            case 8: yAxisLabel = "Mbm"; break;
            case 9: yAxisLabel = "Mby"; break;
            default: yAxisLabel = ""; break;
        }

        graph->GetYaxis()->SetTitle("Observed limit");
        graph->GetXaxis()->SetTitle(yAxisLabel + " (TeV)");

        // Update the canvas
        canvas->Update();

        // Save the plot as a PDF file
        TString pdfFileName = TString::Format("%s.png", yAxisLabel.Data());
        canvas->SaveAs(pdfFileName);

        // Clear the canvas for the next plot
        canvas->Clear();
    }

    // Close the ROOT file
    file.Close();
}

