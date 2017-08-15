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
#include "TCanvas.h"
#include "TLegend.h"



class Canvas
{
protected:
	std::shared_ptr<TCanvas> canv;
	std::shared_ptr<TMultiGraph> mg;
public:
	std::vector<TGraph*> Graphs;
	
	Canvas(const std::string name, const std::string title, int height = 500, int width = 500)
	{
		canv = std::make_shared<TCanvas>(name.c_str(), title.c_str(), height, width);
		mg = std::make_shared<TMultiGraph>(); 
	}
	
	void AddGraph(TGraph* g, const std::string title, const std::string drawOptions)
	{
		g->SetTitle(title.c_str());
		Graphs.push_back(g);
		mg->Add(g, drawOptions.c_str());
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
