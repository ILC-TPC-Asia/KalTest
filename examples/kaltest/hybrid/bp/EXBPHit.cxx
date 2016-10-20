#include "EXBPHit.h"
#include "EXBPMeasLayer.h"
#include "TMath.h"

#include <iostream>
#include <iomanip>

using namespace std;

ClassImp(EXBPHit)

//_____________________________________________________________________
//  ----------------------------------
//  Ctors and Dtor
//  ----------------------------------

EXBPHit::EXBPHit(Int_t m)
        : TVTrackHit(m)
{
}
      
EXBPHit::EXBPHit(const EXBPMeasLayer &ms,
                       Double_t      *x,
                       Double_t      *dx, 
                 const TVector3      &xx,
                       Double_t       b,
                       Int_t          m)
        : TVTrackHit(ms, x, dx, b, m),
          fXX(xx)
{
}

EXBPHit::~EXBPHit()
{
}

//_____________________________________________________________________
//  ----------------------------------
//  Implementation of public methods  
//  ----------------------------------

TKalMatrix EXBPHit::XvToMv(const TVector3 &xv, Double_t /*t0*/) const
{
   const EXBPMeasLayer &ms
                 = dynamic_cast<const EXBPMeasLayer &>(GetMeasLayer());
   TKalMatrix h  = ms.XvToMv(xv);
   Double_t r    = ms.GetR();
   Double_t phih = h(0,0) / r;
   Double_t phim = (*this)(0,0) / r;
   Double_t dphi = phih - phim;

   static Double_t kPi    = TMath::Pi();
   static Double_t kTwoPi = 2 * kPi;

   while (dphi < -kPi) dphi += kTwoPi;
   while (dphi >  kPi) dphi -= kTwoPi;

   h(0,0)  = r * (phim + dphi);
   h(1,0) += 0.;

   return h;
}

void EXBPHit::DebugPrint(Option_t *) const
{
   cerr << "------------------- Site Info -------------------------" << endl;

   for (Int_t i=0; i<GetDimension(); i++) {
      Double_t x  = (*this)(i,0);
      Double_t dx = (*this)(i,1);
      cerr << " x[" << i << "] = " << setw(8) << setprecision(5) << x
           << "    "
           << "dx[" << i << "] = " << setw(6) << setprecision(2) << dx
           << setprecision(7)
           << resetiosflags(ios::showpoint)
           << endl;
   }
   cerr << "-------------------------------------------------------" << endl;
}
