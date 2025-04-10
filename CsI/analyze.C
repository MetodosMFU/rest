
TH1F *spectra = nullptr;
TH1F *xProjection = nullptr;
TH1F *yProjection = nullptr;
TH1F *zProjection = nullptr;
TH2F *xzProjection = nullptr;
TH2F *yzProjection = nullptr;
TH1F *dEdX = nullptr;

std::map <std::string, TH1F*> processHisto;
std::map <std::string, TH1F*> processEnergyHisto;

void analyze (const std::string &fileName){

  TRestRun* run = new TRestRun(fileName);

  run->PrintMetadata();

  TRestGeant4Metadata* G4Metadata = static_cast<TRestGeant4Metadata*>(run->GetMetadataClass("TRestGeant4Metadata"));
  const auto sensitiveVolumeName = G4Metadata->GetSensitiveVolume();
  const auto sensitiveVolID = G4Metadata->GetGeant4GeometryInfo().GetIDFromVolume(sensitiveVolumeName);

  TRestGeant4Event* g4Event = static_cast<TRestGeant4Event*>(run->GetInputEvent());
  double maxX=0, maxY=0, maxZ=0;
  double minX=100000, minY=100000, minZ=100000;
  double EMaxDep = 0;
  const double precision = 0.1;//mm
  for(int i=0;i<run->GetEntries();i++){
      run->GetEntry(i);
      const double sensitiveVolumeEnergy = g4Event->GetEnergyInVolume(sensitiveVolumeName.Data());

      if(sensitiveVolumeEnergy > EMaxDep)EMaxDep = sensitiveVolumeEnergy;
      const auto hits = g4Event->GetHitsInVolume(sensitiveVolID);
      for (size_t hitIndex = 0; hitIndex < hits.GetNumberOfHits(); hitIndex++) {
        auto pos = hits.GetPosition(hitIndex);
        if(pos.X() < minX ) minX= pos.X();
        if(pos.Y() < minY ) minY= pos.Y();
        if(pos.Z() < minZ ) minZ= pos.Z();
        if(pos.X() > maxX ) maxX = pos.X();
        if(pos.Y() > maxY ) maxY = pos.Y();
        if(pos.Z() > maxZ ) maxZ = pos.Z();
      }
  }

  const double maxdX = sqrt((maxX-minX)*(maxX-minX)+ (maxY-minY)*(maxY-minY) +(maxZ-minZ)*(maxZ-minZ) );
  //cout<<"Max "<<maxX<<" "<<maxY<<" "<<maxZ<<endl;
  //cout<<"Min "<<minX<<" "<<minY<<" "<<minZ<<endl;
  //cout<<"MaxDx "<<maxdX <<endl;

  spectra = new TH1F("Spectra","Spectra",1000,0,EMaxDep*1.1 );
  int nBinsX = (maxX-minX)*1.1/precision;
  xProjection = new TH1F("xProjection","xProjection",nBinsX, minX - abs(minX)*0.1, maxX*1.1 );
  int nBinsY = (maxZ-minZ)*1.1/precision;
  yProjection = new TH1F("yProjection","yProjection",nBinsY, minY - abs(minY)*0.1, maxY*1.1 );
  int nBinsZ = (maxZ-minZ)*1.1/precision;
  zProjection = new TH1F("zProjection","zProjection",nBinsX, minZ - abs(minZ)*0.1, maxZ*1.1 );
  xzProjection = new TH2F("xzProjection","xzProjection",nBinsZ, minZ - abs(minZ)*0.1, maxZ*1.1,nBinsX,minX - abs(minX)*0.1, maxX*1.1 );
  yzProjection = new TH2F("yzProjection","yzProjection",nBinsZ,minZ - abs(minZ)*0.1, maxZ*1.1,nBinsY,minY - abs(minY)*0.1, maxY*1.1 );
  int nBins = maxdX*1.1*2.;
  dEdX = new TH1F("dEdX","dEdX",nBins,0, maxdX*1.1 );

  const double nan = std::numeric_limits<double>::infinity();
  std::vector <std::map <std::string, std::pair<int,double> > > processMapVector;
  std::map <std::string, std::pair<int,double> > processMaxMap;

    for(int i=0;i<run->GetEntries();i++){
      run->GetEntry(i);
      std::map <std::string, std::pair<int, double> > processMap;

      const double sensitiveVolumeEnergy = g4Event->GetEnergyInVolume(sensitiveVolumeName.Data());

      spectra->Fill(sensitiveVolumeEnergy);

      auto track = g4Event->GetTrackByID(12);
      if(track == nullptr)continue;
      
      
      const auto hits = track->GetHits();
      double initEn = track->GetInitialKineticEnergy();
      TVector3 prevPos = {0,0,0};
      double dZ=0;

        for (size_t hitIndex = 0; hitIndex < hits.GetNumberOfHits(); hitIndex++) {

          const std::string process = hits.GetProcessName(hitIndex).Data();

          if(process == "Init" || process == "Transportation")continue;

          const TVector3 currentPos = hits.GetPosition(hitIndex);
          xProjection->Fill(currentPos.X());
          yProjection->Fill(currentPos.Y());
          zProjection->Fill(currentPos.Z());
          xzProjection->Fill(currentPos.Z(),currentPos.X());
          yzProjection->Fill(currentPos.Z(),currentPos.Y());
          const double energy = initEn - hits.GetKineticEnergy(hitIndex);
          dZ += (currentPos-prevPos).Mag();
          dEdX->Fill(dZ, energy);

          initEn -= energy;
          prevPos = currentPos;

            auto pM = processMap.find(process);
              if(pM == processMap.end()){
                processMap[process] = std::make_pair(1, energy);
              } else{
                auto p = pM->second; 
                processMap[process] = std::make_pair(p.first+1, p.second+energy);
              }          
        }

          if(!processMap.empty()){
            processMapVector.push_back(processMap);
               for(const auto& [process, p ] : processMap){
                 auto pMax = processMaxMap.find(process);
                 if(pMax == processMaxMap.end() ){
                   processMaxMap[process] = p;
                 } else {
                   auto [maxN, maxEn] = pMax->second;
                   auto [n, energy] = p;
                   bool update = false;
                     if (energy > maxEn ){
                       maxEn = energy;
                       update = true;
                     }
                     if(n > maxN){
                       maxN =  n;
                       update = true;
                     }
                   if(update){
                     //cout<<process<<" "<<maxN<<" "<<maxEn<<endl;
                     processMaxMap[process] = std::make_pair(maxN, maxEn);
                   }
                 }
               }
          }
    }

    for (const auto& processMap : processMapVector){
      for(const auto& [process, p ] : processMap){
        if(processHisto.find(process) == processHisto.end() ){
          auto [maxN, maxEn ] = processMaxMap[process];
          cout<<process<<" "<<maxN << " "<<maxEn <<endl;
          processHisto[process] = new TH1F (process.c_str(), process.c_str(), maxN+3, 0, maxN+2);
          std::string enName = process + "En";
          processEnergyHisto[process] = new TH1F (enName.c_str(), enName.c_str(),1000,0,maxEn*1.1);
        }
        auto [n, energy] = p;
        processHisto[process]->Fill(n);
        processEnergyHisto[process]->Fill(energy);
      }
    }

  TCanvas *spcCan = new TCanvas("spcCan","spcCan");
  spcCan->cd();
  spectra->GetXaxis()->SetTitle("Energy (keV)");
  spectra->GetYaxis()->SetTitle("Counts");
  spectra->Draw();
  TCanvas *proj1DCan = new TCanvas("proj1DCan","proj1DCan");
  proj1DCan->Divide(3,1);
  proj1DCan->cd(1);
  xProjection->GetXaxis()->SetTitle("X (mm)");
  xProjection->Draw();
  proj1DCan->cd(2);
  yProjection->GetXaxis()->SetTitle("Y (mm)");
  yProjection->Draw();
  proj1DCan->cd(3);
  zProjection->GetXaxis()->SetTitle("Z (mm)");
  zProjection->Draw();
  TCanvas *proj2DCan = new TCanvas("proj2DCan","proj2DCan");
  proj2DCan->Divide(2,1);
  proj2DCan->cd(1);
  xzProjection->GetXaxis()->SetTitle("Z (mm)");
  xzProjection->GetYaxis()->SetTitle("X (mm)");
  xzProjection->Draw();
  proj2DCan->cd(2);
  yzProjection->GetXaxis()->SetTitle("Z (mm)");
  yzProjection->GetYaxis()->SetTitle("Y (mm)");
  yzProjection->Draw();

  const int nHistos = processHisto.size();

  TCanvas *dEdXCan = new TCanvas("dEdXCan","dEdXCan");
  dEdX->Scale(1./spectra->GetEntries());
  dEdX->GetXaxis()->SetTitle("Distance (mm)");
  dEdX->GetYaxis()->SetTitle("dE/dX (keV/mm)");
  dEdX->Draw();

  TCanvas *processCan = new TCanvas("processCan","processCan");
  processCan->Divide(nHistos,1);

  int c =1;
  for(const auto& [process, histo] : processHisto){
    processCan->cd(c);
    histo->GetXaxis()->SetTitle("Number of processes per event");
    histo->Draw();
    c++;
  }

  TCanvas *processEnCan = new TCanvas("processEnCan","processEnCan");
  processEnCan->Divide(processEnergyHisto.size(),1);

  c=1;
  for(const auto& [process, histo] : processEnergyHisto){
    processEnCan->cd(c);
    histo->GetXaxis()->SetTitle("Energy per process per event (keV)");
    histo->Draw();
    c++;
  }

}

void saveHisto(TH1 *histo, const std::string &fileName){

  ofstream fq (fileName, ios::out);

  for (int b = 1; b<=histo->GetNbinsX(); b++){
    fq << histo->GetBinCenter(b)<< " "<< histo->GetBinContent(b) << "\n";
  }

  fq.close();

}

