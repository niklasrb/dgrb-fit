#ifndef CANVASWRAPPER_H
#define CANVASWRAPPER_H


#include <memory>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

#include "TROOT.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"


/// A little class that saves the pointer to the canvas
/// and to the graphs on it
class Canvas
{
protected:
	std::shared_ptr<TCanvas> canv;		// shared pointer for easy memory management
	std::shared_ptr<TMultiGraph> mg;
public:
	
	Canvas(const std::string name, const std::string title, int height = 500, int width = 500)
	{
		canv = std::make_shared<TCanvas>(name.c_str(), title.c_str(), height, width);
		mg = std::make_shared<TMultiGraph>(); 
		//Draw();
	}
	
	void AddGraph(TGraph* g, const std::string title, const std::string drawOptions)
	{
		g->SetTitle(title.c_str());
		mg->Add(g, drawOptions.c_str());
		Draw();
	}
	
	void SetxLimits(double xmin, double xmax)
	{
		mg->GetHistogram()->GetXaxis()->SetRangeUser(xmin, xmax);	
		mg->GetXaxis()->SetLimits(xmin, xmax);
		canv->Modified();
		canv->Update();
	}
	
	void SetyLimits(double ymin, double ymax)
	{
		mg->GetHistogram()->GetYaxis()->SetRangeUser(ymin, ymax);
		mg->GetYaxis()->SetLimits(ymin, ymax);
		canv->Modified();
		canv->Update();
	}
	
	void SetAxesTitle(std::string xAxis, std::string yAxis)
	{
		mg->GetXaxis()->SetTitle(xAxis.c_str());
		mg->GetYaxis()->SetTitle(yAxis.c_str());
		canv->Modified();
		canv->Update();
	}
	
	void Draw(std::string drawOptions = "A")
	{
		canv->cd();
		mg->Draw(drawOptions.c_str());
	}
	
	TCanvas& operator()()
	{
		return *canv.get();
	}
	
	TMultiGraph& MultiGraph()
	{
		return *mg.get();
	}
};



#endif
