int read_trigger_definitions(int runNumber = 13078009)
{
  // Load libraries
  // base libraries
  //gSystem->Load("St_base");
  //gSystem->Load("StChain");
  //gSystem->Load("StUtilities");
  //gSystem->Load("StIOMaker");
  //gSystem->Load("StarClassLibrary");
  gROOT->Macro("loadMuDst.C");
  gROOT->Macro("LoadLogger.C");


  // db-related libraries
  gSystem->Load("St_Tables");
  gSystem->Load("StDbLib");
  gSystem->Load("StDbBroker");
  gSystem->Load("St_db_Maker");

  gSystem->Load("StEmcRawMaker");
  gSystem->Load("StEmcAdctoEMaker");
  gSystem->Load("StEEmcUtil");
  gSystem->Load("StEEmcDbMaker");

  gSystem->Load("StTriggerUtilities");

  ifstream in(Form("beginTimes/%d.beginTimes.txt", runNumber));
  if(!in) {
    cout<<"can't open input file\n";
    return 0;
  }
  char run[16];
  char date[16];
  char time[16];
  in>>date>>time;
  cout<<"beginTime:"<<date<<" "<<time<<endl;
  int yr, mnth, dy;
  int hh,mm,ss;
  sscanf(date,"%d-%d-%d", &yr, &mnth, &dy);
  sscanf(time,"%d:%d:%d", &hh, &mm, &ss);

  int idate = yr*10000 + mnth*100 + dy;
  int itime = hh*10000 + mm *100 + ss;
  cout<<idate<<" "<<itime<<endl;

  St_db_Maker *dbMk=new St_db_Maker("db", "MySQL:StarDb", "$STAR/StarDb");
  dbMk->SetDebug();
  dbMk->SetDateTime(idate, itime);
  //  dbMk->SetDateTime(20120730,1); // event or run start time, set to your liking
  dbMk->SetFlavor("ofl");

  dbMk->Init();
  dbMk->Make();

  TDataSet *DB = 0;
  DB = dbMk->GetDataBase("Calibrations/trg/triggerDefinition");
  if (!DB) {
    std::cout << "ERROR: no table found in db, or malformed local db config" << std::endl;
  }

  St_triggerDefinition *dataset = 0;
  dataset = (St_triggerDefinition*) DB->Find("triggerDefinition");
  Int_t rows = dataset->GetNRows();
  if (rows > 1) {
    std::cout << "INFO: found INDEXED table with " << rows << " rows" << std::endl;
  }

  if (dataset) {
    TDatime val[2];
    dbMk->GetValidity((TTable*)dataset,val);
    std::cout << "Dataset validity range: [ " << val[0].GetDate() << "." << val[0].GetTime() << " - " 
	      << val[1].GetDate() << "." << val[1].GetTime() << " ] "
	      << std::endl;

    triggerDefinition_st *table = dataset->GetTable();
    for (Int_t i = 0; i < rows; i++) {
      // sample output of first member variable
      std::cout << i << "th row : " << table[i].comments << std::endl;
      TBufferFile buf(TBuffer::kRead);
      buf.SetBuffer(table[i].trigdef, sizeof(table[i].trigdef), false);
      TObjArray *a = 0;
      buf >> a;
      a->Print();

      cout << "Total entries: " << a->GetEntriesFast() << endl;

      for (int i = 0; i < a->GetEntriesFast(); ++i) {
	StTriggerDefinition* trgdef = (StTriggerDefinition*)a->At(i);
	trgdef->print();
      }
      
    }
  } else {
    std::cout << "ERROR: dataset does not contain requested table" << std::endl;
  }

  return 1;
}
