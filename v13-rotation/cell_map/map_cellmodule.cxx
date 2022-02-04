#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include <any>
#include <map>
#include <string>
#include "TGraph.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TList.h"
#include "TStyle.h"
#include "TText.h"

using namespace std;

void map_cellmodule(){
  double moduler_ = 85.;
  double gap_ = 1.5;
  double nCellRHeight_ = 35.;
  double cellR_ = 2 * moduler_ / nCellRHeight_;
  double C_PI = 3.14159265358979323846;
  
  Double_t a = 2*85./35.;  // cell long radius (or hexagon flat side length)
  Double_t xCenter = 0.;
  Double_t yCenter = 0.;
  Double_t numberOfCellsFlatSide = 12;
  Double_t x[6], y[6];
  Double_t xLoop, yLoop, xTemp, yTemp;
  int yColNum;
  
  map<int, pair<double, double>> cellModulePositionMap_;
  
  for (unsigned moduleID = 0; moduleID < 7; moduleID++) {
    if (moduleID == 0) {
      //============================================================================
      double moduleX = 0.;
      double moduleY = 0.;
      //----------------------------------------------------------------------------
      int ecalMapID = 0;
      xLoop = xCenter - a * 18.0;
      yLoop = yCenter + a * 18.0 / TMath::Sqrt(3);
      for (int xCounter = 0; xCounter < numberOfCellsFlatSide;
           xCounter++) {
        yTemp = yLoop;
        yColNum = numberOfCellsFlatSide + xCounter + 1;
        for (int yCounter = 0; yCounter < yColNum; yCounter++) {
          x[0] = xLoop;
          y[0] = yTemp;
          x[1] = x[0] + a / 2.0;
          y[1] = y[0] + a * TMath::Sqrt(3) / 2.0;
          x[2] = x[1] + a;
          y[2] = y[1];
          x[3] = x[2] + a / 2.0;
          y[3] = y[0];
          x[4] = x[2];
          y[4] = y[3] - a * TMath::Sqrt(3) / 2.0;
          x[5] = x[1];
          y[5] = y[4];
           
          double cellx = (x[0] + x[3]) / 2;
          double celly = (y[1] + y[4]) / 2;
          int cellModuleID = 10 * ecalMapID + moduleID;
          cellModulePositionMap_[cellModuleID] = pair<double, double>(moduleX + cellx, moduleY + celly);
          ecalMapID++;
          yTemp -= a * TMath::Sqrt(3);
        }
        yLoop += a * TMath::Sqrt(3) / 2.0;
        xLoop += 1.5 * a;
      }
      
      // Section 2: Right trapezoid (12cell short base, 23cell long base, 12cell
      xLoop = xCenter;
      yLoop = yCenter + a * 11.0 * TMath::Sqrt(3);
      for (int xCounter = 0; xCounter < numberOfCellsFlatSide; xCounter++) {
        yTemp = yLoop;
        yColNum = 2 * numberOfCellsFlatSide - xCounter - 1;
        for (int yCounter = 0; yCounter < yColNum; yCounter++) {
          x[0] = xLoop;
          y[0] = yTemp;
          x[1] = x[0] + a / 2.0;
          y[1] = y[0] + a * TMath::Sqrt(3) / 2.0;
          x[2] = x[1] + a;
          y[2] = y[1];
          x[3] = x[2] + a / 2.0;
          y[3] = y[0];
          x[4] = x[2];
          y[4] = y[3] - a * TMath::Sqrt(3) / 2.0;
          x[5] = x[1];
          y[5] = y[4];
           
          double cellx = (x[0] + x[3]) / 2;
          double celly = (y[1] + y[4]) / 2;
          int cellModuleID = 10 * ecalMapID + moduleID;
          cellModulePositionMap_[cellModuleID] = pair<double, double>(moduleX + cellx, moduleY + celly);
          ecalMapID++;
          yTemp -= a * TMath::Sqrt(3);
        }
        yLoop -= a * TMath::Sqrt(3) / 2.0;
        xLoop += 1.5 * a;
      }
      //===========================================================================
    }
    
    else {
      //============================================================================
      double moduleX = (2. * moduler_ + gap_) * cos((moduleID-1) * (C_PI / 3.));
      double moduleY = - (2. * moduler_ + gap_) * sin((moduleID-1) * (C_PI / 3.));
      //-----------------------------------------------------------------------------
      int ecalMapID = 0;
      xLoop = xCenter - a * 18.0;
      yLoop = yCenter + a * 18.0 / TMath::Sqrt(3);
      for (int xCounter = 0; xCounter < numberOfCellsFlatSide;
           xCounter++) {
        yTemp = yLoop;
        yColNum = numberOfCellsFlatSide + xCounter + 1;
        for (int yCounter = 0; yCounter < yColNum; yCounter++) {
          x[0] = xLoop;
          y[0] = yTemp;
          x[1] = x[0] + a / 2.0;
          y[1] = y[0] + a * TMath::Sqrt(3) / 2.0;
          x[2] = x[1] + a;
          y[2] = y[1];
          x[3] = x[2] + a / 2.0;
          y[3] = y[0];
          x[4] = x[2];
          y[4] = y[3] - a * TMath::Sqrt(3) / 2.0;
          x[5] = x[1];
          y[5] = y[4];
           
          double cellx = (x[0] + x[3]) / 2;
          double celly = (y[1] + y[4]) / 2;
          int cellModuleID = 10 * ecalMapID + moduleID;
          cellModulePositionMap_[cellModuleID] = pair<double, double>(moduleX + cellx, moduleY + celly);
          ecalMapID++;
          yTemp -= a * TMath::Sqrt(3);
        }
        yLoop += a * TMath::Sqrt(3) / 2.0;
        xLoop += 1.5 * a;
      }
      
      // Section 2: Right trapezoid (12cell short base, 23cell long base, 12cell
      xLoop = xCenter;
      yLoop = yCenter + a * 11.0 * TMath::Sqrt(3);
      for (int xCounter = 0; xCounter < numberOfCellsFlatSide; xCounter++) {
        yTemp = yLoop;
        yColNum = 2 * numberOfCellsFlatSide - xCounter - 1;
        for (int yCounter = 0; yCounter < yColNum; yCounter++) {
          x[0] = xLoop;
          y[0] = yTemp;
          x[1] = x[0] + a / 2.0;
          y[1] = y[0] + a * TMath::Sqrt(3) / 2.0;
          x[2] = x[1] + a;
          y[2] = y[1];
          x[3] = x[2] + a / 2.0;
          y[3] = y[0];
          x[4] = x[2];
          y[4] = y[3] - a * TMath::Sqrt(3) / 2.0;
          x[5] = x[1];
          y[5] = y[4];
           
          double cellx = (x[0] + x[3]) / 2;
          double celly = (y[1] + y[4]) / 2;
          int cellModuleID = 10 * ecalMapID + moduleID;
          cellModulePositionMap_[cellModuleID] = pair<double, double>(moduleX + cellx, moduleY + celly);
          ecalMapID++;
          yTemp -= a * TMath::Sqrt(3);
        }
        yLoop -= a * TMath::Sqrt(3) / 2.0;
        xLoop += 1.5 * a;
      }
      //==============================================================================
    }
  }
  
  ofstream file;
  file.open ("cellmodule.txt");
  for (auto const& [cellModuleID, cellCenter] : cellModulePositionMap_) {
    file << cellModuleID << " " << cellCenter.first << " " << cellCenter.second << "\n";
  }
  file.close();
  
  const Int_t n_points = cellModulePositionMap_.size();
  Double_t x_pos[n_points], y_pos[n_points];
  
  int ID_pos = 0;
  for (auto const& [cellModuleID, cellCenter] : cellModulePositionMap_) {
    x_pos[ID_pos] = cellCenter.first;
    y_pos[ID_pos] = cellCenter.second;
    ID_pos++;
  }
  
  // Create the canvas and graph
  TCanvas* mycanvas = new TCanvas("mycanvas","Ecal Flower face",900,900);
  mycanvas->cd();
  TGraph *myGraph = new TGraph(n_points, x_pos, y_pos);
  myGraph->SetTitle("Ecal Flower face Cell Center Positions");
  myGraph->SetMarkerStyle(4);
  myGraph->SetMarkerColor(kBlack);
  myGraph->Draw("AP"); // A for axes, P for polyline, E for errors
  mycanvas->SaveAs("cellModulePositions.pdf");
  
  return;
}

