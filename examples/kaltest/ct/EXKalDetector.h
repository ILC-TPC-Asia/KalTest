#ifndef EXDETECTOR_H
#define EXDETECTOR_H

#include "TVector3.h"         // from ROOT
#include "TVKalDetector.h"    // from KalTrackLib
#include "EXMeasLayer.h"

class EXKalDetector : public TVKalDetector {
public:
   // Ctor and Dtor
   EXKalDetector(Int_t m = 100);
   ~EXKalDetector();

   // Utility methods

   static Double_t GetBfield (const TVector3 &xx = TVector3())
                               { return fgBfield; }

private:
   static Double_t fgBfield;   // magnetic field [kG]

   ClassDef(EXKalDetector,1)   // Sample hit class
};

#endif
