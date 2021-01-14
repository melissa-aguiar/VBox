void Muons::UserInit(Long64_t nentries)
{
    h0=new TH1F*[128];
    h1=new TH1F*[128];
    h2=new TH1F*[128];
    h3=new TH1F*[128];
    c = new TCanvas*[128];
    c2 = new TCanvas*[128];

    for (int side=0;side<2;side++){
        for(int mod=0;mod<64;mod++){
            sprintf(buf0, "Hist_lado%d/mod%d",side,mod);
            h0[64*side+mod]=new TH1F(buf0,buf0,100,-2000,8000);
            sprintf(buf1, "Hist_lado%d/mod%d",side,mod);
            h1[64*side+mod]=new TH1F(buf1,buf1,100,-2000,8000);
            sprintf(buf2, "Hist_lado%d/mod%d",side,mod);
            h2[64*side+mod]=new TH1F(buf2,buf2,100,-2000,8000);
            sprintf(buf3, "Hist_lado%d/mod%d",side,mod);
            h3[64*side+mod]=new TH1F(buf3,buf3,100,-2000,8000);
            for(int v=0;v<1000;v++){
                ROC_pds[side][mod][v]=0;
                ROC_fas[side][mod][v]=0;
                ROC_pdd[side][mod][v]=0;
                ROC_fad[side][mod][v]=0;
            }
        }
    }
}

void Muons::leitura_pesos()
{
    ifstream inCoefText;
   	
    inCoefText.open("Pesos_filtro_ajust.txt");  
   
    Double_t vCoef[11];
    int side, mod,ch; 
  
    if(!inCoefText.is_open()){ exit(EXIT_FAILURE); }
	while (1)    
	{
		inCoefText >> vCoef[0] >> vCoef[1] >> vCoef[2] >> vCoef[3] >> vCoef[4] >> vCoef[5] >> vCoef[6]>> vCoef[7]>> vCoef[8] >> vCoef[9] >> vCoef[10];
		if (!inCoefText.good()) break;
			side=vCoef[0];
			mod=vCoef[1];
			ch=vCoef[2];
            for (unsigned int j=0;j<7;j++){
                peso[side][mod][ch][j] = vCoef[j+4];		
                //cout<<peso[side][mod][ch][j]<<" ";
			}
		    //cout<<endl;
	}
	inCoefText.close();	
}

void Muons::leitura_calibracao()
{
    ifstream coefcalib;

    coefcalib.open("calibration_2017_2_ajust.txt");
    int lado, modulo, canal;
    if(!coefcalib.is_open()){exit(EXIT_FAILURE); }
    while (1)
    {
        for (int j=0;j<5;j++){
            coefcalib >> caltxt[j];
        }
        
        if (!coefcalib.good()) break;
            lado=caltxt[0];
            modulo=caltxt[1];
            canal=caltxt[2];

            for (int i=0; i<2; i++){
                cal[lado][modulo][canal][i]=caltxt[i+3];
                //cout << cal[lado][modulo][canal][i] << " ";
            }
            //cout << endl;
    }
    coefcalib.close();
}

void Muons::leitura_coef_den()
{
    ifstream coefdenoising;

    coefdenoising.open("Coeficientes_Denoising.txt");
    int l, m;
    Double_t coefs_den[6];
    if(!coefdenoising.is_open()){exit(EXIT_FAILURE); }
    while (1)
    {
        for (int j=0;j<6;j++){
            coefdenoising >> coefs_den[j];
        }
        
        if (!coefdenoising.good()) break;
            l=coefs_den[0];
            m=coefs_den[1];

            for (int i=0; i<4; i++){
                coefs[l][m][i]=coefs_den[i+2];
                //cout << coefs[l][m][i] << " ";
            }
            //cout << endl;
    }
    coefdenoising.close();
}

void Muons::Est_MFMuon()
{
    
    for(int side=0;side<2;side++){
        for(int mod=0;mod<64;mod++){
            for(int ch=0;ch<4;ch++){
                mf_muon[side][mod][ch] = 0;
                for(int i=0;i<7;i++){
                    mf_muon[side][mod][ch] += SmpMuon[side][mod][ch][i]*peso[side][mod][ch][i]; 
                    //cout<<peso[side][mod][ch][i]<<" ";      
                }
                //cout << endl;
                //cout << mf[side][mod][ch] << " ";
                mev_muon[side][mod][ch] = (mf_muon[side][mod][ch]*cal[side][mod][ch][0])+cal[side][mod][ch][1];
                /*for (int i=0; i<2; i++){
                    cout << cal[side][mod][ch][i] << " ";
                }
                cout << endl;*/
                //cout << mev[side][mod][ch] << " ";
            }
            est_s_muon[side][mod][0] = mev_muon[side][mod][0]+mev_muon[side][mod][1];
            est_s_muon[side][mod][1] = mev_muon[side][mod][2]+mev_muon[side][mod][3];
            //cout << est_soma[side][mod][0] << " " << est_soma[side][mod][1];
            //cout << endl;
        }
    }
}

void Muons::Est_MFRuido()
{
    
    for(int side=0;side<2;side++){
        for(int mod=0;mod<64;mod++){
            for(int ch=0;ch<4;ch++){
                mf_ruido[side][mod][ch] = 0;
                for(int i=0;i<7;i++){
                    mf_ruido[side][mod][ch] += SmpNoise[side][mod][ch][i]*peso[side][mod][ch][i]; 
                    //cout<<peso[side][mod][ch][i]<<" ";      
                }
                //cout << endl;
                //cout << mf[side][mod][ch] << " ";
                mev_ruido[side][mod][ch] = (mf_ruido[side][mod][ch]*cal[side][mod][ch][0])+cal[side][mod][ch][1];
                /*for (int i=0; i<2; i++){
                    cout << cal[side][mod][ch][i] << " ";
                }
                cout << endl;*/
                //cout << mev[side][mod][ch] << " ";
            }
            est_s_ruido[side][mod][0] = mev_ruido[side][mod][0]+mev_ruido[side][mod][1];
            est_s_ruido[side][mod][1] = mev_ruido[side][mod][2]+mev_ruido[side][mod][3];
            //cout << est_soma[side][mod][0] << " " << est_soma[side][mod][1];
            //cout << endl;
        }
    }
}

void Muons::Est_Denoising_Muon()
{
    for(int mod=0;mod<64;mod++){
        for(int ch=0;ch<4;ch++){
            est_d_muon[0][mod][0] = coefs[0][mod][0]*mev_muon[0][mod][0]+coefs[0][mod][1]*mev_muon[0][mod][1];
            est_d_muon[0][mod][1] = coefs[0][mod][2]*mev_muon[0][mod][2]+coefs[0][mod][3]*mev_muon[0][mod][3];
            est_d_muon[1][mod][0] = coefs[1][mod][0]*mev_muon[1][mod][0]+coefs[1][mod][1]*mev_muon[1][mod][1];
            est_d_muon[1][mod][1] = coefs[1][mod][2]*mev_muon[1][mod][2]+coefs[1][mod][3]*mev_muon[1][mod][3];
            /*cout << est_d_muon[0][mod][0] << " " << est_d_muon[0][mod][1];
            cout << endl;*/
        }
    }
}

void Muons::Est_Denoising_Ruido()
{
    for(int mod=0;mod<64;mod++){
        for(int ch=0;ch<4;ch++){
            est_d_ruido[0][mod][0] = coefs[0][mod][0]*mev_ruido[0][mod][0]+coefs[0][mod][1]*mev_ruido[0][mod][1];
            est_d_ruido[0][mod][1] = coefs[0][mod][2]*mev_ruido[0][mod][2]+coefs[0][mod][3]*mev_ruido[0][mod][3];
            est_d_ruido[1][mod][0] = coefs[1][mod][0]*mev_ruido[1][mod][0]+coefs[1][mod][1]*mev_ruido[1][mod][1];
            est_d_ruido[1][mod][1] = coefs[1][mod][2]*mev_ruido[1][mod][2]+coefs[1][mod][3]*mev_ruido[1][mod][3];
                /*cout << est_den[0][mod][0] << " " << est_den[0][mod][1];
                cout << endl;*/
        }
    }
}

void Muons::Eficiencia()
{
    TMatrixD rms_soma(1, 2);
    TMatrixD rms_den(1, 2);
    for (unsigned int par=0;par<2;par++){
        for(int mod=0;mod<64;mod++){
            rms_soma.Zero();
            rms_den.Zero();
            /*for (int i=0;i<2;i++){
                    rms_soma[0][i] = coefs[par][mod][1]*/
        }
    }
}

void Muons::Fill_hist()
{
    for (int side=0;side<2;side++){
        for(int mod=0;mod<64;mod++){
            h0[side*64+mod]->Fill(est_s_ruido[side][mod][0]+est_s_ruido[side][mod][1]);
            h1[side*64+mod]->Fill(est_d_ruido[side][mod][0]+est_d_ruido[side][mod][1]);
            h2[side*64+mod]->Fill(est_s_muon[side][mod][0]+est_s_muon[side][mod][1]);
            h3[side*64+mod]->Fill(est_d_muon[side][mod][0]+est_d_muon[side][mod][1]);
        }
    }
}

void Muons::save_hist()
{
    for (int side=0;side<2;side++){
        for(int mod=0;mod<64;mod++){
            sprintf(buf4, "Canvas_lado%d/mod%d",side,mod);
            c[64*side+mod]=new TCanvas[buf4,"Hist_distr", 1];
            c[side*64+mod]->cd(1);
            h0[side*64+mod]->Draw();
            h0[side*64+mod]->GetXaxis()->SetTitle("Estimacao [MeV]");
            h0[side*64+mod]->SetLineColor(kRed);
            h1[side*64+mod]->Draw("Same");
            h1[side*64+mod]->SetLineColor(kBlue);
            h2[side*64+mod]->Draw("Same");
            h2[side*64+mod]->SetLineColor(kGreen);
            h3[side*64+mod]->Draw("Same");
            h3[side*64+mod]->SetLineColor(kBlack);
            auto* legend = new TLegend(0.70 , 0.70, .95, .95);
            media0=h0[side*64+mod]->GetMean();
            variancia0=h0[side*64+mod]->GetRMS();
            sprintf(leg0, "Ruido sem denoising, media=%d  RMS= %d", media0, variancia0);
            media1=h1[side*64+mod]->GetMean();
            variancia1=h1[side*64+mod]->GetRMS();
            sprintf(leg1, "Ruido com denoising, media=%d  RMS= %d", media1, variancia1);
            media2=h2[side*64+mod]->GetMean();
            variancia2=h2[side*64+mod]->GetRMS();
            sprintf(leg2, "Sinal sem denoising, media=%d  RMS= %d", media2, variancia2);
            media3=h3[side*64+mod]->GetMean();
            variancia3=h3[side*64+mod]->GetRMS();
            sprintf(leg3, "Sinal com denoising, media=%d  RMS= %d", media3, variancia3);
            legend->AddEntry(h0[side*64+mod],leg0);
            legend->AddEntry(h1[side*64+mod],leg1);
            legend->AddEntry(h2[side*64+mod],leg2);
            legend->AddEntry(h3[side*64+mod],leg3);
            legend->Draw();
            sprintf(filename,"GRAFICO_side_%d_mod_%d_.pdf",side,mod);
            c[side*64+mod]->Print(filename);
            c[side*64+mod]->Close();
        }
    }
}

void Muons::detec_falsoa()
{    
    for (int side=0;side<2;side++){
        for(int mod=0;mod<64;mod++){
            ruido[side][mod]=est_s_ruido[side][mod][0]+est_s_ruido[side][mod][1];
            r_den[side][mod]=est_d_ruido[side][mod][0]+est_d_ruido[side][mod][1];
            muon[side][mod]=est_s_muon[side][mod][0]+est_s_muon[side][mod][1];
            m_den[side][mod]=est_d_muon[side][mod][0]+est_d_muon[side][mod][1];
            for(int v=0;v<1000;v++){
                if(muon[side][mod]>(-1000+6*v)){ROC_pds[side][mod][v]=ROC_pds[side][mod][v]+1;}
                if(ruido[side][mod]>(-1000+6*v)){ROC_fas[side][mod][v]=ROC_fas[side][mod][v]+1;}
                if(m_den[side][mod]>(-1000+6*v)){ROC_pdd[side][mod][v]=ROC_pdd[side][mod][v]+1;}
                if(r_den[side][mod]>(-1000+6*v)){ROC_fad[side][mod][v]=ROC_fad[side][mod][v]+1;}
            }         
        }
    }
}

void Muons::plot_ROC()
{
    for (int side=0;side<2;side++){
        for(int mod=0;mod<64;mod++){
            for(int v=0;v<1000;v++){
                ROC_pds[side][mod][v]=ROC_pds[side][mod][v]/nentries;
                ROC_fas[side][mod][v]=ROC_fas[side][mod][v]/nentries;
                ROC_pdd[side][mod][v]=ROC_pdd[side][mod][v]/nentries;
                ROC_fad[side][mod][v]=ROC_fad[side][mod][v]/nentries;
            }
        }
    }

    Int_t n=1000;
    for (int s=0;s<2;s++){
        for(int m=0;m<64;m++){
            TGraph *roc[2];
            TMultiGraph *mroc = new TMultiGraph();
            sprintf(buf5, "Canvas_lado%d/mod%d",s,m);
            c2[64*s+m] = new TCanvas[buf5,"Curva ROC", 1];
            
            roc[0] = new  TGraph(n,ROC_fas[s][m],ROC_pds[s][m]); 
            roc[0]->SetName("gr4");
			roc[0]->SetTitle("Sem Denoising");
			roc[0]->SetDrawOption("a");
			roc[0]->SetLineColor(2);
			roc[0]->SetLineWidth(1);
			roc[0]->SetFillStyle(0);
			roc[0]->SetMarkerStyle(20);
			roc[0]->SetMarkerColor(2);
			mroc->Add(roc[0]);
            
            roc[1] = new  TGraph(n,ROC_fad[s][m],ROC_pdd[s][m]);
            roc[1]->SetName("gr3");
			roc[1]->SetTitle("Com Denoising"); //Legenda
			roc[1]->SetDrawOption("a");
			roc[1]->SetLineColor(4); //Cor da linha
			roc[1]->SetLineWidth(1); //Espessura da linha
			roc[1]->SetFillStyle(0);
			roc[1]->SetMarkerStyle(20); //Define a bolinha como marcador
			roc[1]->SetMarkerColor(4); //Cor da bolinha
			mroc->Add(roc[1]);

            mroc->Draw("APL");
            mroc->GetXaxis()->SetTitle("False Alarm (%)");
            mroc->GetYaxis()->SetTitle("Efficiency (%)");
            //mgroc->SetMinimum(60);
            //mgroc->SetMaximum(100);
            mroc->GetYaxis()->SetRangeUser(0.97,0.99);
            mroc->GetXaxis()->SetRangeUser(0,0.2);
            gPad->Modified();
            c2[64*s+m]->BuildLegend(0.9,0.9,0.68,0.78); 
            char filename7[100];
            sprintf(filename7,"ROC_sid_%d_mod_%d.pdf",s,m);
            c2[64*s+m]->Print(filename7);
        }
    }
}

void Muons::SalveMuons()
{
    ofstream coef0;
    coef0.open("Sample_Muons");
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

