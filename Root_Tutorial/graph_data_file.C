void graph_data_file()
{
	TGraph *gr = new TGraph();
	
	gr->SetMarkerStyle(kFullCircle);
	gr->SetLineWidth(2);
	gr->SetLineColor(kRed);
	
	gr->SetTitle("Graph");
	gr->GetXaxis()->SetTitle("X Values");
	gr->GetYaxis()->SetTitle("Y Values");
	
	fstream file;
	file.open("data2.txt", ios::in);
	
	while(1)
	{
		double x, y;
		file >> x >> y;
		gr->SetPoint(gr->GetN(), x, y); //GetN gets the next datapoint
		if(file.eof()) break;
	}
	file.close();
	
	TCanvas *c1 = new TCanvas();
	gr->Draw("ALP"); //L means all the data points are connected with a line, * means the points are stars
}
