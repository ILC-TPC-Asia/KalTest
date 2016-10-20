#ifndef EXITDETECTOR_H
#define EXITDETECTOR_H

#include "EXVKalDetector.h"

class EXITKalDetector : public EXVKalDetector {
public:
   EXITKalDetector(Int_t m = 100);
   ~EXITKalDetector();

   ClassDef(EXITKalDetector,1)   // Sample hit class
};

#endif
