void hist_read_input()
{
	TH1F *hist = new TH1F("hist", "My Histogram", 7, 1, 7);
	
	fstream file;
	file.open("data.txt", ios::in);
	
	double value;
	
	while(1)
	{
		file >> value;
		hist->Fill(value);
		if(file.eof()) break; //if its the end of the file, break
	}
	
	file.close();
	
	hist->GetXaxis()->SetTitle("Grade");
	hist->GetYaxis()->SetTitle("Entries");
	
	TCanvas *c1 = new TCanvas();
	hist->Draw();
}
