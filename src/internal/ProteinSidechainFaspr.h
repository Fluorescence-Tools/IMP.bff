#ifndef IMPBFF_FASPR_H
#define IMPBFF_FASPR_H

// -------- from FasprAAName.h --------
/* Vendored 1:1 from FASPR (github.com/tommyhuangthu/FASPR, version
 * 20200309) into IMP.bff by tools/port_faspr.py -- behaviour-identical port
 * wrapped in namespace IMP::bff::faspr. Only mechanical changes: namespace,
 * include guards/prefixes, exit() -> throwing faspr_fail(), MSVC pragmas
 * dropped, sprintf -> snprintf. See IMP/bff/ProbeRotamerLibrary.h for the license
 * preserved below and the citation. */
/* FASPR original header follows verbatim. */
/*******************************************************************************************************************************
This file is a part of the protein side-chain packing software FASPR

Copyright (c) 2020 Xiaoqiang Huang (tommyhuangthu@foxmail.com, xiaoqiah@umich.edu)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation 
files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, 
modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the 
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
********************************************************************************************************************************/
#include <string>

namespace IMP {
namespace bff {
namespace faspr {


using namespace std;

char Three2One(string aa3);
void One2Three(char aa1, string &aa3);


} // namespace faspr
} // namespace bff
} // namespace IMP

// -------- from FasprUtility.h --------
/* Vendored 1:1 from FASPR (github.com/tommyhuangthu/FASPR, version
 * 20200309) into IMP.bff by tools/port_faspr.py -- behaviour-identical port
 * wrapped in namespace IMP::bff::faspr. Only mechanical changes: namespace,
 * include guards/prefixes, exit() -> throwing faspr_fail(), MSVC pragmas
 * dropped, sprintf -> snprintf. See IMP/bff/ProbeRotamerLibrary.h for the license
 * preserved below and the citation. */
/* FASPR original header follows verbatim. */
/*******************************************************************************************************************************
This file is a part of the protein side-chain packing software FASPR

Copyright (c) 2020 Xiaoqiang Huang (tommyhuangthu@foxmail.com, xiaoqiah@umich.edu)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation 
files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, 
modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the 
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
********************************************************************************************************************************/
#include <cmath>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <stdexcept>

namespace IMP {
namespace bff {
namespace faspr {


using namespace std;

typedef vector<float> FV1;
typedef vector<vector<float> > FV2;
typedef vector<vector<vector<float> > > FV3;
typedef vector<vector<vector<vector<float> > > > FV4;
typedef vector<double> DV1;
typedef vector<vector<double> > DV2;
typedef vector<vector<vector<double> > > DV3;
typedef vector<vector<vector<vector<double> > > > DV4;
typedef vector<int> IV1;
typedef vector<vector<int> > IV2;
typedef vector<vector<vector<int> > > IV3;
typedef vector<vector<vector<vector<int> > > > IV4;
typedef vector<string> SV1;
typedef vector<vector<string> > SV2;

const float DELTA=1.0e-15;
const float BONDDIST=2.1;
const float PI=3.1415926;
const float DEG2RAD=PI/180.0e0;
const float RAD2DEG=180.0e0/PI;

float Distance(FV1 &p1,FV1 &p2);
float Angle(FV1 &p1,FV1 &p2,FV1 &p3);
float Dihedral(FV1 &p1,FV1 &p2,FV1 &p3,FV1 &p4);
bool Internal2Cartesian(FV1 &c1,FV1 &c2,FV1 &c3,FV1 &p,FV1 &cc);


float VectorDotProduct(FV1 &c1,FV1 &c2);
void VectorCrossProduct(FV1 &c1,FV1 &c2,FV1 &cc);
void VectorAdd(FV1 &c1,FV1 &c2,FV1 &cc);
void VectorSubtract(FV1 &c1,FV1 &c2,FV1 &cc);
void VectorMinus(FV1 &cc);
float VectorAngle(FV1 &c1,FV1 &c2);
bool VectorNormalization(FV1 &c);
void VectorMultiply(float par,FV1 &c);
float VectorModulo(FV1 &c);
float Sign(float c);
void MatrixByVector(FV2 &mtr,FV1 &vec,FV1 &cc);
bool VectorL2Norm(FV1 &data);


//! FASPR's exit(0) error paths become exceptions in-process.
[[noreturn]] inline void faspr_fail() {
  throw std::runtime_error(
      "IMP.bff FASPR port: fatal error (see stderr)");
}

} // namespace faspr
} // namespace bff
} // namespace IMP

// -------- from FasprStructure.h --------
/* Vendored 1:1 from FASPR (github.com/tommyhuangthu/FASPR, version
 * 20200309) into IMP.bff by tools/port_faspr.py -- behaviour-identical port
 * wrapped in namespace IMP::bff::faspr. Only mechanical changes: namespace,
 * include guards/prefixes, exit() -> throwing faspr_fail(), MSVC pragmas
 * dropped, sprintf -> snprintf. See IMP/bff/ProbeRotamerLibrary.h for the license
 * preserved below and the citation. */
/* FASPR original header follows verbatim. */
/*******************************************************************************************************************************
This file is a part of the protein side-chain packing software FASPR

Copyright (c) 2020 Xiaoqiang Huang (tommyhuangthu@foxmail.com, xiaoqiah@umich.edu)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation 
files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, 
modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the 
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
********************************************************************************************************************************/
#include <iomanip>
#include <fstream>

namespace IMP {
namespace bff {
namespace faspr {


using namespace std;

struct Residue
{
  string name;
  char chID;
  int pos;
  char ins;
  SV1 atNames;
  string atTypes;
  FV2 xyz;
};

typedef vector<Residue> PV1;
typedef vector<vector<Residue> > PV2;

class Structure
{
public:
  string seq;
  PV1 pdb;
  int nres;

  ~Structure();
  void Pdb2Fas();
  void ReadPDB(string &pdbfile);
  void OutputPDB(PV1 &pdb);
  void OutputPDB(PV1 &pdb,string &pdbfile);
  void WritePDB(string &pdbfile);
};



} // namespace faspr
} // namespace bff
} // namespace IMP

// -------- from FasprRotamerBuilder.h --------
/* Vendored 1:1 from FASPR (github.com/tommyhuangthu/FASPR, version
 * 20200309) into IMP.bff by tools/port_faspr.py -- behaviour-identical port
 * wrapped in namespace IMP::bff::faspr. Only mechanical changes: namespace,
 * include guards/prefixes, exit() -> throwing faspr_fail(), MSVC pragmas
 * dropped, sprintf -> snprintf. See IMP/bff/ProbeRotamerLibrary.h for the license
 * preserved below and the citation. */
/* FASPR original header follows verbatim. */
/*******************************************************************************************************************************
This file is a part of the protein side-chain packing software FASPR

Copyright (c) 2020 Xiaoqiang Huang (tommyhuangthu@foxmail.com, xiaoqiah@umich.edu)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation 
files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, 
modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the 
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
********************************************************************************************************************************/
#include <map>
#include <cstring>

namespace IMP {
namespace bff {
namespace faspr {


using namespace std;

struct Topology
{
  int nchi;
  SV1 atnames;
  IV2 former3AtomIdx;
  FV2 ic;
};

class RotamerBuilder:public Structure
{
public:
  ~RotamerBuilder();
  void LoadSeq();
  void LoadSeq(string &seqfile);
  void LoadParameter();
  void LoadBBdepRotlib2010();
  void BuildSidechain();
  void RotlibFromBinary2Text(string binlibfile,string &txtlibfile);
  void RotlibFromText2Binary(string &fulltextlib,string &binlibfile);
  void PhiPsi();
  void AssignSidechainTopology();
  int LoadBackbone(int site);
  void SideChain(int site,int rot,FV2& rxyz);

  IV1 subStat;
  IV1 nrots;
  FV1 phi,psi;
  FV3 chi;
  FV2 probRot;
  FV1 maxProb;
  PV1 stru;
  FV4 sc;
  map<char,float> wRotlib;
  map<char,Topology> sidechainTopo;
};


} // namespace faspr
} // namespace bff
} // namespace IMP

// -------- from FasprSelfEnergy.h --------
/* Vendored 1:1 from FASPR (github.com/tommyhuangthu/FASPR, version
 * 20200309) into IMP.bff by tools/port_faspr.py -- behaviour-identical port
 * wrapped in namespace IMP::bff::faspr. Only mechanical changes: namespace,
 * include guards/prefixes, exit() -> throwing faspr_fail(), MSVC pragmas
 * dropped, sprintf -> snprintf. See IMP/bff/ProbeRotamerLibrary.h for the license
 * preserved below and the citation. */
/* FASPR original header follows verbatim. */
/*******************************************************************************************************************************
This file is a part of the protein side-chain packing software FASPR

Copyright (c) 2020 Xiaoqiang Huang (tommyhuangthu@foxmail.com, xiaoqiah@umich.edu)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation 
files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, 
modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the 
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
********************************************************************************************************************************/

namespace IMP {
namespace bff {
namespace faspr {


#define WGT_HBOND       1.0
#define WGT_SSBOND      6.0
#define CACB_DIST_CUT   2.35
#define RESI_DIST_CUT   4.25
#define SEC_CUT        15.0
//atomic parameters
#define RADIUS_C1       1.78
#define RADIUS_C2       1.40
#define RADIUS_C3       2.30
#define RADIUS_C4       1.91
#define RADIUS_C5       1.96
#define RADIUS_C6       1.73
#define RADIUS_C7       1.43
#define RADIUS_C8       1.99
#define RADIUS_O1       1.48
#define RADIUS_O2       1.44
#define RADIUS_O3       1.40
#define RADIUS_O4       1.43
#define RADIUS_N1       1.42
#define RADIUS_N2       1.69
#define RADIUS_N3       1.56
#define RADIUS_N4       1.70
#define RADIUS_S1       2.15
#define RADIUS_S2       1.74
#define DEPTH_C1        0.25
#define DEPTH_C2        0.14
#define DEPTH_C3        0.30
#define DEPTH_C4        0.37
#define DEPTH_C5        0.48
#define DEPTH_C6        0.38
#define DEPTH_C7        0.07
#define DEPTH_C8        0.38
#define DEPTH_O1        0.22
#define DEPTH_O2        0.27
#define DEPTH_O3        0.07
#define DEPTH_O4        0.12
#define DEPTH_N1        0.08
#define DEPTH_N2        0.24
#define DEPTH_N3        0.46
#define DEPTH_N4        0.48
#define DEPTH_S1        0.44
#define DEPTH_S2        0.40
//VDW
#define DSTAR_MIN_CUT   0.015
#define DSTAR_MAX_CUT   1.90
#define VDW_REP_CUT    10.0
//Hbond energy
#define OPT_HBOND_DIST  2.8
#define MIN_HBOND_DIST  2.6
#define MAX_HBOND_DIST  3.2
#define MIN_HBOND_THETA 90.
#define MIN_HBOND_PHI   90.
//SSbond energy
#define OPT_SSBOND_DIST 2.03
#define MIN_SSBOND_DIST 1.73
#define MAX_SSBOND_DIST 2.53
#define OPT_SSBOND_ANGL 105.
#define MIN_SSBOND_ANGL 75.
#define MAX_SSBOND_ANGL 135.

using namespace std;

class SelfEnergy:public RotamerBuilder
{
public:
  IV1 bestrot;
  FV2 eTableSelf;
  IV2 conMap;
  FV2 radius,depth;
  IV2 atomIdx;

  ~SelfEnergy();
  void AssignConMap();

  void SetVdwPar();
  float VDWType(int a,int b,float rij,float dist);
  float VDWEnergyAtomAndAtom(float ddash,float eij);
  float RotamerPreferenceEnergy(int site,int rot);

  void EnergySidechainAndBackbone(int site,FV1 &ener);
  void EnergyRotamerSidechainAndFixedSidechain(int site,FV1 &ener);
  void CalcSelfEnergy();
  float HbondEnergyAtomAndAtom(FV1 &DBxyz,FV1 &Dxyz,FV1 &Axyz,FV1 &ABxyz,float Dangle,float Aangle);
  float SSbondEnergyAtomAndAtom(FV1 &CA1xyz,FV1 &CB1xyz,FV1 &SG1xyz,FV1 &SG2xyz,FV1 &CB2xyz,FV1 &CA2xyz);
  float EnergyPolarSidechainAndBackbone(int site1,int rot1,int site2);
  float EnergyPolarSidechainAndSidechain(int site1,int rot1,int site2,int rot2);
};


} // namespace faspr
} // namespace bff
} // namespace IMP

// -------- from FasprPairEnergy.h --------
/* Vendored 1:1 from FASPR (github.com/tommyhuangthu/FASPR, version
 * 20200309) into IMP.bff by tools/port_faspr.py -- behaviour-identical port
 * wrapped in namespace IMP::bff::faspr. Only mechanical changes: namespace,
 * include guards/prefixes, exit() -> throwing faspr_fail(), MSVC pragmas
 * dropped, sprintf -> snprintf. See IMP/bff/ProbeRotamerLibrary.h for the license
 * preserved below and the citation. */
/* FASPR original header follows verbatim. */
/*******************************************************************************************************************************
This file is a part of the protein side-chain packing software FASPR

Copyright (c) 2020 Xiaoqiang Huang (tommyhuangthu@foxmail.com, xiaoqiah@umich.edu)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation 
files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, 
modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the 
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
********************************************************************************************************************************/

namespace IMP {
namespace bff {
namespace faspr {


using namespace std;

void Transpose(float **&tab,int i,int j);

class PairEnergy:public SelfEnergy{
public:
  float**** eTablePair;

  ~PairEnergy();
  void CalcPairEnergy();
  bool EnergyRotamerSidechainAndRotamerSidechain(int site1,int site2,float **&tab);
  void ShowPairEnergy();
  void ShowPairEnergy(int site1, int site2);
};



} // namespace faspr
} // namespace bff
} // namespace IMP

// -------- from FasprSearch.h --------
/* Vendored 1:1 from FASPR (github.com/tommyhuangthu/FASPR, version
 * 20200309) into IMP.bff by tools/port_faspr.py -- behaviour-identical port
 * wrapped in namespace IMP::bff::faspr. Only mechanical changes: namespace,
 * include guards/prefixes, exit() -> throwing faspr_fail(), MSVC pragmas
 * dropped, sprintf -> snprintf. See IMP/bff/ProbeRotamerLibrary.h for the license
 * preserved below and the citation. */
/* FASPR original header follows verbatim. */
/*******************************************************************************************************************************
This file is a part of the protein side-chain packing software FASPR

Copyright (c) 2020 Xiaoqiang Huang (tommyhuangthu@foxmail.com, xiaoqiah@umich.edu)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation 
files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, 
modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the 
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
********************************************************************************************************************************/
#include <set>
#include <stack>
#include <ctime>

namespace IMP {
namespace bff {
namespace faspr {


using namespace std;

#define FASPR_DEE_THRESHOLD 0.0

/**************************************************************************
                     RankElement and RankSize
***************************************************************************/
template <class Elem> class RankElement{
public:
  RankElement(int a,Elem b):idx(a),element(b){}
  int idx;
  Elem element;
  bool operator < (const RankElement &m)const {
    return element<m.element;
  }
  bool operator > (const RankElement &m)const {
    return element>m.element;
  }
};


//sort elements according to the size
template <class Elem> class RankSize{
public:
  RankSize(int a,Elem b):idx(a),element(b){}
  int idx;
  Elem element;
  bool operator < (const RankSize &m)const {
    return element.size()<m.element.size();
  }
  bool operator > (const RankSize &m)const {
    return element.size()>m.element.size();
  }
};

/**************************************************************************
                                 Graph
The whole residue interaction graph consists of several subgraphs/clusters, 
each cluster can be represented as a map:
  map(current Vertex)->neighbor Vertices
***************************************************************************/
typedef map<int, set <int> > Graph;
void ShowGraph(Graph &graph);


/**************************************************************************
                                 Bag
***************************************************************************/
typedef enum{
  Root,
  Inner,
  Leaf,
  None
}BagType;

class Bag{
public:
  set<int> left;
  set<int> right;
  set<int> total;
  int parentBagIdx;
  set<int> childBagIdx;
  int type;
  int childCounter;

  /* data structure to record solution*/
  IV1 lsites;
  IV1 rsites;
  IV1 tsites;
  IV2 lrots;
  IV2 rrots;
  FV2 Etcom;
  IV2 Rtcom;
  IV1 indices;
  bool deployFlag;

  /* for backtrack*/
  float EGMEC;
  IV1 RLGMEC;
  IV1 RRGMEC;

  Bag(){
    type=Leaf;
    deployFlag=false;
  }
  void ShowBag();
};


/**************************************************************************
                         Tree Decomposition
***************************************************************************/
class TreeDecomposition{
public:
  vector <Bag> bags;
  vector <Bag> connBags;

  void Subgraph2TreeDecomposition(int index,Graph &graph);
  void MergeBags(int depth);
  int CheckTreewidth();

};


class Solution:public PairEnergy{
public:
  ~Solution();
  IV1 unfixres;
  bool DEESearch(IV1 &pos);
  int DEEGoldstein(IV1& pos);
  int DEEsplit(IV1& pos);
  void Pick(int site,int rot);

  vector< Graph > graphs;
  void ConstructAdjMatrix(int nunfix,IV2 &adjMatrix);
  void ConstructSubgraphs(int nunfix,IV1 &visited,IV2 &adjMatrix,IV2 &flagMatrix);
  void FindSubgraphByDFS(Graph &graph,int u,IV1 &visited,IV2 &adjMatrix,IV2 &flagMatrix,stack<int> &vertices);
  void ShowGraphs();
  void GraphEdgeDecomposition(IV2 &adjMatrix,float threshold);

  TreeDecomposition tree;
  void TreeDecompositionBottomToTopCalcEnergy();
  void CalcRightBagRotamerCombinationEnergy(Bag &leafbag,int depth,float &Etmp,IV1 &Rtmp,FV1 &Ercom,IV2 &Rrcom);
  void CalcLeftBagRotamerCombinationEnergy(Bag &rootbag,int depth,float &Etmp,IV1 &Rtmp,FV1 &Elcom,IV2 &Rlcom);
  void GetLeftBagRotamerCombination(Bag &leafbag,int depth,IV1 &Rtmp,IV2 &Rlcom);
  void BagDeploySites(Bag &leafbag);
  void LeafBagCalcEnergy(Bag &leafbag,IV2 &Rlcom);
  void CombineChildIntoParentBag(Bag &leafbag,Bag &parbag,IV2 &Rclcom);
  void SubsetCheck(IV1 &subset,IV1 &fullset, IV1 &indices);
  void RootBagFindGMEC(Bag &rootbag);
  void TreeDecompositionTopToBottomAssignRotamer(Bag &parbag,Bag &childbag);
  void TreeDecompositionRelease();
  void Search();
};


} // namespace faspr
} // namespace bff
} // namespace IMP

#endif  // IMPBFF_FASPR_H
