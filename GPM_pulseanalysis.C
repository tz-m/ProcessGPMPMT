void GPM_pulseanalysis(std::string filename)
{
  gROOT->Reset();

  TFile * infile = TFile::Open(filename.c_str(),"READ");
  
  Double_t calibration_an = 1.10848e12;
  Double_t calibration_top = 8.9663e11;
  Double_t elem_chg = 1.602e-19;

  TTreeReader reader("pulsetree",infile);
  TTreeReaderValue<Int_t> channel(reader,"channel");
  TTreeReaderValue<Float_t> amplitude(reader,"amplitudeVolt");

  Double_t minchargetop = 9e99, minchargean = 9e99;
  Double_t maxchargetop = -9e99, maxchargean = -9e99;
  while (reader.Next())
    {
      if (*channel == 7)
	{
	  Double_t q = *amplitude/(calibration_top*elem_chg);
	  if (q > maxchargetop) maxchargetop = q;
	  if (q < minchargetop && q > 1) minchargetop = q;
	}
      else if (*channel == 8)
	{
	  Double_t q = *amplitude/(calibration_an*elem_chg);
	  if (q > maxchargean) maxchargean = q;
	  if (q < minchargean && q > 1) minchargean = q;
	}
    }
  reader.SetEntry(0);

  RooRealVar elec_an("elec_an","Q_{anode}",minchargean,maxchargean);//1.0/(calibration_an*elem_chg));
  RooRealVar elec_top("elec_top","Q_{top}",minchargetop,maxchargetop);//1.0/(calibration_top*elem_chg));

  RooRealVar mean_an("mean_an","mean",0);
  RooRealVar gain_an("gain_an","gain",1e6,1e4,1e7);
  RooRealVar m_an("m_an","m",0.1,0.0001,10);

  RooRealVar mean_top("mean_top","mean",0);
  RooRealVar gain_top("gain_top","gain",1e6,1e4,1e7);
  RooRealVar m_top("m_top","m",0.1,0.0001,10);

  RooGenericPdf polya_an("polya_an","polya_an","pow(m_an,m_an)*(1/ROOT::Math::tgamma(m_an))*(1/gain_an)*pow((elec_an-mean_an)/gain_an,m_an-1)*exp(-m_an*(elec_an-mean_an)/gain_an)",RooArgSet(elec_an,mean_an,gain_an,m_an));

  RooGenericPdf polya_top("polya_top","polya_top","pow(m_top,m_top)*(1/ROOT::Math::tgamma(m_top))*(1/gain_top)*pow((elec_top-mean_top)/gain_top,m_top-1)*exp(-m_top*(elec_top-mean_top)/gain_top)",RooArgSet(elec_top,mean_top,gain_top,m_top));

  RooDataSet data_an("data_an","data_an",elec_an);
  RooDataSet data_top("data_top","data_top",elec_top);

  TH1F * anodes = new TH1F("anodes","",256,0,0);
  TH1F * tops = new TH1F("tops","",256,0,0);

  while (reader.Next())
    {
      if (*channel == 7)
	{
	  elec_top = *amplitude / (calibration_top*elem_chg);
	  data_top.add(elec_top);
	  tops->Fill(elec_top.getValV());
	}
      else if (*channel == 8)
	{
	  elec_an = *amplitude / (calibration_an*elem_chg);
	  data_an.add(elec_an);
	  anodes->Fill(elec_an.getValV());
	}
    }

  Float_t maxtop = tops->GetBinContent(tops->GetMaximumBin());
  elec_top.setRange("fitrange_top",tops->GetBinCenter(tops->FindLastBinAbove(0.9*maxtop)),tops->GetBinCenter(tops->FindLastBinAbove(0.05*maxtop)));
  Float_t maxan = anodes->GetBinContent(anodes->GetMaximumBin());
  elec_an.setRange("fitrange_an",anodes->GetBinCenter(anodes->FindLastBinAbove(0.9*maxan)),anodes->GetBinCenter(anodes->FindLastBinAbove(0.05*maxan)));

  
  RooFitResult* r_an = polya_an.fitTo(data_an,RooFit::Save(kTRUE),RooFit::Range("fitrange_an"));
  RooFitResult* r_top = polya_top.fitTo(data_top,RooFit::Save(kTRUE),RooFit::Range("fitrange_top"));

  RooPlot* frame_an = elec_an.frame(RooFit::Title("Anode Charge"));
  data_an.plotOn(frame_an,RooFit::DataError(RooAbsData::SumW2),RooFit::Binning(256));
  polya_an.plotOn(frame_an);
  TPaveLabel * tan = new TPaveLabel(0.6,0.8,0.95,0.88,Form("gain = %.2e +/- %.2e",gain_an.getVal(),gain_an.getError()),"NDC");
  tan->SetFillColor(0);
  tan->SetBorderSize(1);
  frame_an->addObject(tan);
  frame_an->getAttText()->SetTextSize(0.3);

  RooPlot* frame_top = elec_top.frame(RooFit::Title("Top Charge"));
  data_top.plotOn(frame_top,RooFit::DataError(RooAbsData::SumW2),RooFit::Binning(256));
  polya_top.plotOn(frame_top);
  TPaveLabel * ttop = new TPaveLabel(0.6,0.8,0.95,0.88,Form("gain = %.2e +/- %.2e",gain_top.getVal(),gain_top.getError()),"NDC");
  ttop->SetFillColor(0);
  ttop->SetBorderSize(1);
  frame_top->addObject(ttop);
  frame_top->getAttText()->SetTextSize(0.3);

  std::cout << "!!!!!!!!!!!!!!!!!!! ANODES !!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  r_an->Print();
  std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" << std::endl;
  std::cout << "!!!!!!!!!!!!!!!!!!! TOPS !!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  r_top->Print();
  std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" << std::endl;

  TCanvas* canv = new TCanvas("canv","GainView",2400,1200);
  canv->Divide(2);
  
  canv->cd(1); 
  gPad->SetLeftMargin(0.15); 
  frame_an->GetYaxis()->SetTitleOffset(1.8); 
  frame_an->Draw();
  
  canv->cd(2); 
  gPad->SetLeftMargin(0.15); 
  frame_top->GetYaxis()->SetTitleOffset(1.8); 
  frame_top->Draw();

  canv->SaveAs("gain.png");

}

