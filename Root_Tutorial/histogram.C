void histogram()
{
	TH1F *hist = new TH1F("hist", "My Histogram", 100, 0, 100);
	
	hist->Fill(10); //fill it with the value 10
	hist->Fill(90);
	
	hist->GetXaxis()->SetTitle("X Axis");
	hist->GetYaxis()->SetTitle("Y Axis");
	
	TCanvas *c1 = new TCanvas();
	hist->Draw();
}
