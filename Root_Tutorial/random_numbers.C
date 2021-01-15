void random_numbers()
{
	 TRandom2 *rand = new TRandom2(1); //if its empty, the seed is set automatically. if 0, the seed is random
	 
	 TH1F *hist = new TH1F("hist", "Histogram", 100, 0, 100);
	 
	 for(int i=0; i<10000; i++)
	 {
	 	double r = rand->Rndm()*100;
	 	cout << r << endl;
	 	hist->Fill(r);
	 }
	 
	 TCanvas *c1 = new TCanvas();
	 hist->GetYaxis()->SetRangeUser(0, 500);
	 hist->Draw();
}
