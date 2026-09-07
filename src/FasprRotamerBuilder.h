/* Vendored 1:1 from FASPR (github.com/tommyhuangthu/FASPR, version
 * 20200309) into IMP.bff by tools/port_faspr.py -- behaviour-identical port
 * wrapped in namespace IMP::bff::faspr. Only mechanical changes: namespace,
 * include guards/prefixes, exit() -> throwing faspr_fail(), MSVC pragmas
 * dropped, sprintf -> snprintf. See IMP/bff/RotamerLibrary.h for the license
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
#ifndef FASPR_ROTAMER_BUILDER_H
#define FASPR_ROTAMER_BUILDER_H

#include "FasprStructure.h"
#include <fstream>
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
#endif
