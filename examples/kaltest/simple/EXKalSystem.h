#ifndef EXKALSYSTEM_H
#define EXKALSYSTEM_H
//*************************************************************************
//* ===================
//*  EXKalSystem Class
//* ===================
//*
//* (Description)
//*   Sample steering class for Kalman filter
//* (Requires)
//*     TVKalSystem
//* (Provides)
//*     class EXKalSystem
//* (Update Recored)
//*   2003/09/30  Y.Nakashima	Original version.
//*
//*************************************************************************
                                                                                
#include "TVKalSystem.h"

//_____________________________________________________________________
//  ------------------------------
//  Sample Kalman Filter class
//  ------------------------------
//
                                                                                
class EXKalSystem : public TVKalSystem {
public:
   EXKalSystem(Int_t n = 1);
   ~EXKalSystem() {} 

private:

   ClassDef(EXKalSystem,1)  // Base class for Kalman Filter
};
                                                                                
//=======================================================
// inline functions
//=======================================================
                                                                                
#endif
