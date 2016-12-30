#ifndef _CalcClashScore_include_
#define _CalcClashScore_include_

# include "../../LibraryFiles/TypeDef.h"
# include "../../LibraryFiles/minimize.h"
# include "XRayStruct.h"
# include "XRay.h"
# include "LatticeStruct.h"

void calcClashScoreRecord(vector<AtomStruct> &Atoms1, vector<AtomStruct> &Atoms2, Real &clashCC, Real &clashCS, Real &clashSS);
void calcClashScoreRecord(ProteinStruct &Protein, vector<ProteinStruct> &allProteins, int index, Real cutoff, Real &clashCC, Real &clashCS, Real &clashSS);
Real calcClashScoreRecord(ProteinStruct &Protein, XRayParamStruct &params, Real &clashCC, Real &clashCS, Real &clashSS);
void calcClashScore(AtomStruct &Atom, vector<AtomStruct> &Atoms, Real &clashCC, Real &clashCS, Real &clashSS);
void calcClashScore(vector<AtomStruct> &Atoms1, vector<AtomStruct> &Atoms2, Real &clashCC, Real &clashCS, Real &clashSS);
void calcClashScore(ProteinStruct &Protein, vector<ProteinStruct> &allProteins, int index, Real cutoff, Real &clashCC, Real &clashCS, Real &clashSS);
void calcClashScoreByAtom(AtomStruct &Atom, Vector &center, vector<ProteinStruct> &allProteins, int index, Real cutoff, Real &clashCC, Real &clashCS, Real &clashSS);
void calcClashScoreByAtom(ProteinStruct &Protein, vector<ProteinStruct> &allProteins, int index, Real cutoff, Vector &clashCC, Vector &clashCS, Vector &clashSS);
Real calcClashScoreByAtom(ProteinStruct Protein, XRayParamStruct &params, Vector &clashCC, Vector &clashCS, Vector &clashSS, Vector &clashScoreAtom);
Real calcClashScoreByAtom(ProteinStruct Protein, XRayParamStruct &params, Vector &clashScoreAtom);
Real calcClashScore(ProteinStruct Protein, XRayParamStruct &params, Real &clashCC, Real &clashCS, Real &clashSS);
Real calcClashScore(ProteinStruct &Protein, XRayParamStruct &params);
void printClashScore(ProteinStruct &Protein, XRayParamStruct &params);
Real calcClashScore(vector<ProteinStruct> &Proteins, vector<ProteinStruct> &allProteins, Real &clashCC, Real &clashCS, Real &clashSS, XRayParamStruct &params);
Real calcClashScore(Vector &degOfFreedom, ProteinStruct Protein, XRayParamStruct &params, Real &clashCC, Real &clashCS, Real &clashSS);
Real calcClashScore(Vector &degOfFreedom, ProteinStruct &Protein, XRayParamStruct &params);
void makeUnitCellAndImagesFast(vector<ProteinStruct> &Proteins);
Real calcClashScoreFast(Vector &degOfFreedom, vector<ProteinStruct> &Proteins, XRayParamStruct &params, Real &clashCC, Real &clashCS, Real &clashSS);
Real calcClashScoreFast(Vector &degOfFreedom, vector<ProteinStruct> &Proteins, XRayParamStruct &params);
void printClashScore(Vector &degOfFreedom, ProteinStruct &Protein, XRayParamStruct &params);
Real calcClashScore(vector<ProteinStruct> Proteins, XRayParamStruct &params);
Real calcClashScore(Vector &degOfFreedoms, vector<ProteinStruct> Proteins, XRayParamStruct &params);
#endif
