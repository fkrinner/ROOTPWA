//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      mass dependent fit likelihood rank 1!
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef MASSDEPFITLIKELI_HH
#define MASSDEPFITLIKELI_HH

// Base Class Headers ----------------
#include "Math/IFunction.h"

// Collaborating Class Headers -------
#include <vector>
#include <map>
#include <string>
#include "pwacomponent.h"

// Collaborating Class Declarations --
class TTree;
class TF1;

namespace rpwa {

  class fitResult;



  class massDepFitLikeli : public ROOT::Math::IBaseFunctionMultiDim {
  public:
    
    // Constructors/Destructors ---------
    massDepFitLikeli(){}
    virtual ~massDepFitLikeli(){}
    
    
    // Accessors -----------------------
    virtual unsigned int NDim() const;
    unsigned int NDataPoints() const; /// number of data points in fit

    // Modifiers -----------------------
    void init(TTree* fitresulttree,
	      TF1* finalStatePhaseSpace,
	      pwacompset* compset,
	      double mmin=0, double mmax=5000);
    


    // Operations ----------------------
    virtual double DoEval  (const double* par) const;
    
    virtual IBaseFunctionMultiDim* Clone()const {return new massDepFitLikeli(*this);}
    
    
    
  private:
    
    // Private Data Members ------------
    TTree* _tree;
    TF1* _finalStatePS;
    pwacompset* _compset;
    fitResult* _rhom; // measured spindensity matrix;
    std::vector<std::string> _wlist;
    std::vector<unsigned int> _index; // wave indices
    
    double  _mmin, _mmax; // fitrange
    // Private Methods -----------------
    
  };
  
} // end namespace

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
