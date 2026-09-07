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
#ifndef FASPR_AANAME_H
#define FASPR_AANAME_H

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
#endif
