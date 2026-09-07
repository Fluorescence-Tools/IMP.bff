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
#ifndef FASPR_UTILITY_H
#define FASPR_UTILITY_H

#include <cmath>
#include <iostream>
#include <vector>
#include <string>
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
#endif
