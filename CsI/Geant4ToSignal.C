
constexpr double photoElectronsPerkeV = 3;//Average number of photoelectrons per keV
constexpr double timeResolution = 0.04;//us
constexpr double timeOffset = 4;//us
constexpr int nTimeBins = 512;//Number of timeBins in the readout electronics
constexpr int shapingTime = 10;//units are timeBins
constexpr double gain = 10;
constexpr int Nr = 5 * shapingTime;
constexpr double noiseLevel = 2.;//ADCs
constexpr double ADCOffset = 250;//ADCs

std::vector<double> initShaper(){
  std::vector<double> response(Nr);
    for (int i = 0; i < Nr; i++) {
      const double coeff = ((Double_t)i) / shapingTime;
      response[i] = TMath::Exp(-3. * coeff) * coeff * coeff * coeff * sin(coeff);
    }
  double sum=0;
  for (int n = 0; n < Nr; n++) sum += response[n];
  for (int n = 0; n < Nr; n++) response[n] = response[n] * gain / sum;

  return response;
}

static std::vector<double> responseShaper = initShaper();

void Geant4ToSignal(const std::string &inputFile){

  TRestRun inputRun(inputFile);
  inputRun.PrintMetadata();

  TRestGeant4Metadata* G4Metadata = static_cast<TRestGeant4Metadata*>(inputRun.GetMetadataClass("TRestGeant4Metadata"));
  const auto sensitiveVolumeName = G4Metadata->GetSensitiveVolume();
  const auto sensitiveVolID = G4Metadata->GetGeant4GeometryInfo().GetIDFromVolume(sensitiveVolumeName);

  TRestGeant4Event* g4Event = static_cast<TRestGeant4Event*>(inputRun.GetInputEvent());
  
  TRestRawSignalEvent rawSignalEvent;
  TRestDetectorHitsEvent hitsEvent;

  TRestRun outRun (inputFile);
  outRun.SetRunNumber(inputRun.GetRunNumber());
  std::string runTag = inputRun.GetRunTag().Data();
  outRun.SetRunTag(runTag);
  outRun.SetRunType("RAWSIGNAL");
  outRun.SetOutputFileName("Run[fRunNumber]_[fRunType]_[fRunTag]_[fRunUser].root");
  outRun.FormOutputFile();
  outRun.AddEventBranch(g4Event);
  outRun.AddEventBranch(&rawSignalEvent);
  outRun.SetStartTimeStamp(inputRun.GetStartTimestamp());

  const int pulseStart = std::round(timeOffset/timeResolution);
  TVector2 baselineRange (0,pulseStart);
  TVector2 integralRange (pulseStart,nTimeBins-1);

  const int runEntries = inputRun.GetEntries();
  for(int i=0;i<runEntries;i++){
    inputRun.GetEntry(i);

    double firstHitTime = std::numeric_limits<double>::max();

    const auto nTracks = g4Event->GetNumberOfTracks();
    if(nTracks<=0)continue;
    for (unsigned int trackIndex = 0; trackIndex < nTracks; trackIndex++) {
      const auto trackG4 = g4Event->GetTrackPointer(trackIndex);

      const auto hits = trackG4->GetHits();
      const auto nHits = hits.GetNumberOfHits();
      //cout<<particleName<<" "<<nHits<<endl;
      if(nHits<=0)continue;
        for (unsigned int n = 0; n < nHits; n++) {
          const double energy = hits.GetEnergy(n);
          const double hitTime = hits.GetTime(n);
          if( hits.GetVolumeId(n) != sensitiveVolID)continue;
          if (energy <= 0)continue;
          if(hitTime < firstHitTime)firstHitTime = hitTime;
        }
      }

    int triggerOffset =timeOffset/timeResolution;

    double totalPhotons =0;
    std::map<int, std::vector<int> > rawSignalMap;
    double totalEDepG4 =0;
    double thresholdIntegral = 0;
    double amplitude = 0;
    double averageTime = 0;
    double risetime = 0;
    double baseline = 0;
    double baselineSigma=0;
    double width =0;

    std::vector<int> data (nTimeBins, 0);

    for (unsigned int trackIndex = 0; trackIndex < nTracks; trackIndex++) {
      const auto trackG4 = g4Event->GetTrackPointer(trackIndex);
      const std::string particleName = trackG4->GetParticleName().Data();

      const auto hits = trackG4->GetHits();
      const auto nHits = hits.GetNumberOfHits();
        for (unsigned int n = 0; n < nHits; n++) {
          if( hits.GetVolumeId(n) != sensitiveVolID)continue;
          const double energy = hits.GetEnergy(n);
          if (energy <= 0)continue;
          const TVector3& position = hits.GetPosition(n);
          const double hitTime = hits.GetTime(n);
          const int nPhotons = std::round(energy*photoElectronsPerkeV);
          if(nPhotons<=1)continue;
          const double diffTime = hitTime - firstHitTime;
          int timeBin = std::round(diffTime/timeResolution + triggerOffset);
             if(timeBin< 0 || timeBin >= nTimeBins){
               std::cout<<"Warning time bin is "<<timeBin<<" but maximum bin is "<<nTimeBins<<" "<< firstHitTime<<" "<<hitTime<<" " << diffTime <<endl;
               continue;
             }
          data[timeBin] += nPhotons;
          totalEDepG4 += energy;
          totalPhotons += nPhotons;
         }
        }

        if(totalPhotons<=0)continue;

        double res = 1./TMath::Sqrt(totalPhotons);
        TRandom3 smear(0);
        const double smearFactor = smear.Gaus(1.0, res);
        for (auto& val : data)val *=smearFactor;
        std::vector<int> shapedData(nTimeBins,0);
          for (int m = 0; m < nTimeBins; m++) {
            if(data[m]==0)continue;
              for (int n = 0; m + n < nTimeBins && n < Nr; n++){
                shapedData[m + n] += responseShaper[n] * data[m];
              }
          }

        rawSignalEvent.Initialize();
        rawSignalEvent.SetID(g4Event->GetID());
        rawSignalEvent.SetSubID(g4Event->GetSubID());
        rawSignalEvent.SetTimeStamp(g4Event->GetTimeStamp());
        rawSignalEvent.SetSubEventTag(g4Event->GetSubEventTag());
        TRestRawSignal rawSignal;
        rawSignal.SetSignalID(1);

        TRandom3 noise(0);
            for(const auto& val : shapedData){
              int ADC = val+noise.Gaus(0, noiseLevel) + ADCOffset;
              if(ADC<0 || ADC >std::numeric_limits<short>::max()) ADC =std::numeric_limits<short>::max();
              rawSignal.AddPoint((Short_t)ADC);
            }

        rawSignalEvent.AddSignal(rawSignal);

        rawSignalEvent.SetBaseLineRange(baselineRange);
        rawSignalEvent.SetRange(integralRange);

        const int nSignals = rawSignalEvent.GetNumberOfSignals();
        if(nSignals<=0)continue;

          for (int s = 0; s < nSignals; s++) {
            TRestRawSignal* sgnl = rawSignalEvent.GetSignal(s);
            sgnl->InitializePointsOverThreshold(TVector2(2., 1.),2);
            double integral = sgnl->GetThresholdIntegral();
            risetime += sgnl->GetRiseTime();
            width += sgnl->GetMaxPeakWidth();
            baseline += sgnl->GetBaseLine();
            baselineSigma += sgnl->GetBaseLineSigma();
            double amp = sgnl->GetMaxValue();
            thresholdIntegral += integral;
            amplitude += amp;
            double time = sgnl->GetMaxPeakBin()*timeResolution;
          }

        outRun.GetAnalysisTree()->SetObservableValue("thresholdIntegral", thresholdIntegral);
        outRun.GetAnalysisTree()->SetObservableValue("totalAmplitude", amplitude);
        outRun.GetAnalysisTree()->SetObservableValue("risetime", risetime);
        outRun.GetAnalysisTree()->SetObservableValue("width", width);
        outRun.GetAnalysisTree()->SetObservableValue("baseline", baseline);
        outRun.GetAnalysisTree()->SetObservableValue("baselineSigma", baselineSigma);
        outRun.GetAnalysisTree()->SetObservableValue("totalEDepG4",totalEDepG4);

        outRun.GetAnalysisTree()->SetEventInfo(&rawSignalEvent);
        outRun.GetEventTree()->Fill();
        outRun.GetAnalysisTree()->Fill();

        if(i%1000==0)std::cout<<"Processed "<<i<<" events out of "<< runEntries<<endl;
  }

  outRun.SetEndTimeStamp(inputRun.GetEndTimestamp());
  outRun.UpdateOutputFile();
  G4Metadata->Write(G4Metadata->GetName());
  outRun.CloseFile();
}


