#ifndef CANVASWRAPPER_H
#define CANVASWRAPPER_H

#include "TROOT.h"
#include "TGraph.h"
#include "TCanvas.h"

#include <memory>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>


/// GraphWrapper
/// Combines TGraph and the actual location of the x and y coordinates
/// To ease memory management 
/*
class GraphWrapper
{
protected:
	double* x;
	double* y;
	unsigned int n;
	TGraph graph;
	std::string drawOptions;
	
public:
	GraphWrapper(unsigned int n, double* _x, double* _y, std::string drawOptions);
	GraphWrapper(const GraphWrapper& gw);
	GraphWrapper& operator =(const GraphWrapper& gw);
	~GraphWrapper();
	
	TGraph& operator ->();
	TGraph& Graph();
	void Draw();
};

GraphWrapper::GraphWrapper(unsigned int n, double* _x, double* _y, std::string drawOptions) : n(n) , drawOptions(drawOptions)
{
	x = new double[n];
	y = new double[n];
	for(unsigned int i = 0; i < n; i++)
	{
		x[i] = _x[i];
		y[i] = _y[i];
	}
	graph = TGraph(n, x, y);
}

GraphWrapper::GraphWrapper(const GraphWrapper& gw) : GraphWrapper(gw.n, gw.x, gw.y, gw.drawOptions)
{}

GraphWrapper& GraphWrapper::operator =(const GraphWrapper& gw)
{
	if(&gw == this) return *this;
	GraphWrapper tmp(gw);
	std::swap(tmp.n, n); std::swap(tmp.drawOptions, drawOptions);
	std::swap(tmp.x, x); std::swap(tmp.y, y);
	std::swap(tmp.graph, graph);
	return *this;
}

TGraph& GraphWrapper::operator ->()
{
	return graph;
}
TGraph& GraphWrapper::Graph()
{
	return graph;
}
void GraphWrapper::Draw()
{
	std::cout << drawOptions << '\t' << x << '\t' << y << '\t' << y[2] << std::endl;
	graph.Draw(drawOptions.c_str());
}

GraphWrapper::~GraphWrapper()
{
	delete []x;
	delete []y;
}*/

/// CanvasWrapper
/// Combines Canvas and its Graphs
/// To ease memory management
class CanvasWrapper
{
protected:
	std::shared_ptr<TCanvas> canv;
	
public:
	typedef std::pair<std::shared_ptr<TGraph>, std::string> Graph;
	std::vector<Graph > graphs;
	
	
	CanvasWrapper(const std::string& name, const std::string& title, int width, int height);
	
	void AddGraph(const TGraph& g, const std::string drawOptions);
	void Redraw();
	
	TCanvas& Canvas();
	TCanvas& operator ->();
	
	void SaveToFile(std::string file);
	
};

CanvasWrapper::CanvasWrapper(const std::string& name, const std::string& title, int width = 500, int height = 500)
{
	canv = std::make_shared<TCanvas>(name.c_str(), title.c_str(), width, height);
}
	
void CanvasWrapper::AddGraph(const TGraph& g, const std::string drawOptions)
{
	graphs.push_back(std::make_pair(std::make_shared<TGraph>(g), drawOptions));
	Redraw();
}

void CanvasWrapper::Redraw()
{
	//canv->Clear();
	canv->cd();
	for(auto it = graphs.begin(); it != graphs.end(); ++it)
	{
		std::cout << it->first << '\t' << it->second << std::endl;
		it->first->Draw((it->second + "same" + (it == ++graphs.begin() ? "A" : "")).c_str());
	}
}
	
void CanvasWrapper::SaveToFile(std::string file)
{
	
}

TCanvas& CanvasWrapper::Canvas()
{
	return *canv.get();
}



#endif
