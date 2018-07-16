{

TChain* chain = new TChain("T");

char filename[100];  
FILE* fd;
for (int t=1868;t<=2214;t++) {
  sprintf(filename,"/w/work3504/e01015/peter/offlana/elastics/nd_threshold/e01015phy_12C_4pass_%d.root",t);
  fd = fopen(filename,"r");
  if (fd != NULL) {
        cout << "Including run " << t << endl;
    chain->Add(filename);
    fclose(fd);
      } else {
        cout << "Run " << t << " not here." << endl;
  }
}  

}
