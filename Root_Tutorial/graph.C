void graph()
{
	double x[5] = {1, 2, 3, 4, 5};
	double y[5] = {1, 4, 9, 16, 25};
	
	TGraph *gr = new TGraph(5, x, y);
	
	//gr->SetMarkerStyle(4); //opcional, o 4 faz pontos como bolinhas
	//gr->SetMarkerSize(1); //opcional
	
	TCanvas *c1 = new TCanvas();
	gr->Draw("AL*"); //L means all the data points are connected with a line, * means the points are stars
}
