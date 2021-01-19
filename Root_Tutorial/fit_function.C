void fit_function()
{
	TH1F *hist = new TH1F("hist", "My Gaussian", 100, 0, 10);
	
	TRandom2 *rand = new TRandom2(3);
	
	fstream file;
	file.open("data3.txt", ios::out);
	
	for(int i=0; i<1000; i++)
	{
		double r = rand->Gaus(5,1); //(mean, sigma)
		file << r << endl;
	}
	
	file.close();
	
	file.open("data3.txt", ios::in);
		
	double value;
	
	while(1)
	{
		file >> value;
		hist->Fill(value);
		if(file.eof()) break; //if its the end of the file, break
	}
	
	file.close();
	
	hist->GetXaxis()->SetTitle("Distribution");
	hist->GetYaxis()->SetTitle("Entries");
	
	TF1 *fit = new TF1("fit", "gaus", 0, 10);
	
	fit->SetParameter(0, 40);
	fit->SetParameter(1, 5);
	fit->SetParameter(2, 1);
	
	TCanvas *c1 = new TCanvas();
	hist->Draw();
	
	hist->Fit("fit", "R");
	
	double mean = fit->GetParameter(1);
	double sigma = fit->GetParameter(2);
	
	cout << mean/sigma << endl;
}
