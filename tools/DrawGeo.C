{
	
# include <stdio.h>
# include <iostream.h>
# include <fstream.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
  cout << 1 << endl;
	TGeoManager * gm;
	cout << 2 << endl;
	//Read top volume
	gm = TGeoManager::Import("geo.gdml");
	cout << 3 << endl;
        gm->SetVisLevel(10);
	
	cout << 4 << endl;
			
	TGeoVolume * gvol=gm->GetTopVolume();
	TCanvas *c1 = new TCanvas("c1","c1",100,100,500,500);
	c1->cd();
			
	gvol->Draw();
	c1->Modified();
	c1->Update();
   
	
	
}
