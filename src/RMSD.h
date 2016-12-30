#ifndef _RMSD_include_
#define _RMSD_include_
void moveProteinFractional(ProteinStruct &Protein, Real fracX, Real fracY, Real fracZ);
void correctForArbitraryAxis(ProteinStruct &Protein1, ProteinStruct &Protein2);
Real calcRMSDSymmetryOperation(ProteinStruct Protein1, ProteinStruct &Protein2, int index);
Real calcRMSDReference(ProteinStruct &Protein1, ProteinStruct &Protein2, Real x, Real y, Real z, int &index);
Real calcRMSDSymmetryOperations(ProteinStruct Protein1, ProteinStruct &Protein2, int &bestIndex);
Real calcRMSD(ProteinStruct Protein1, ProteinStruct &Protein2, Real fracX, Real fracY, Real fracZ, int &bestIndex);
void getAlternateOrigins(string spaceGroup, Matrix &origin);
Real calcRMSD(ProteinStruct &Protein1, ProteinStruct &Protein2, string superimposeType);
Real calcRMSD_NoAlign(ProteinStruct &Protein1, ProteinStruct &Protein2);
Real calcRMSD(vector<ProteinStruct> &Proteins1, vector<ProteinStruct> &Proteins2, vector<int> &permutation, string superimposeType);
Real calcPointXtalRmsd(ProteinStruct &Protein1, ProteinStruct &Protein2);
Real calcRMSD(Matrix &rmsdPair, vector<int> &permutation);
void getMatchingAtomsNoChain(vector<AtomStruct> &Atoms1, vector<AtomStruct> &Atoms2, vector<AtomStruct> &AtomsOut);
Real calcRMSD(vector<ProteinStruct> &Proteins1, vector<ProteinStruct> &Proteins2, string superimposeType);
#endif
