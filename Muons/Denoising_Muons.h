void Muons::SalveMuons()
{
    ofstream coef0;
    coef0.open("Sample_Muons.txt");
    //coef0.open("Estimativa_denoising");
    for(int side=0;side<2;side++){
        for(int mod=0;mod<64;mod++){
            coef0 << side << " " << mod << " ";
            for(int ch=0;ch<4;ch++){
                //coef0 << mf[side][mod][ch] << " ";
                for(int i=0;i<7;i++){
                    coef0 << SmpMuon[side][mod][ch][i] << " ";
                }
            }
            coef0 << endl;
        }
    }
    coef0.close();
}

