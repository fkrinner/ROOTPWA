///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
#ifdef __CINT__


#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;


#pragma link C++ class TFitBin+;
#pragma link C++ class rpwa::fitResult+;
#pragma link C++ class rpwa::pwaPlotter+;
#pragma link C++ class rpwa::TPwaFitGraphErrors+;
#pragma link C++ class TCMatrix+;
#pragma link C++ class TCovEllipse+;
#pragma link C++ class TMCMCMeta+;



// TFitResult produces a name clash for ROOT versions from 5.25.0 on
#include "RVersion.h"
// rootcint has problems with this: #if ROOT_VERSION_CODE < ROOT_VERSION(5,25,0)
#if ROOT_VERSION_CODE < 334080  // make sure ROOT version is below 5.25.0
#pragma link C++ class TFitResult+;
#endif 


#endif
