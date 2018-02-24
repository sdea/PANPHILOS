#ifndef SOLVESBMPROCESSORCONTACTANGLE3D_H_INCLUDED
#define SOLVESBMPROCESSORCONTACTANGLE3D_H_INCLUDED

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace plb;

template<typename T>
class solveSBMContactAngle3D : public BoxProcessingFunctional3D {

public:

    // Constructor
    solveSBMContactAngle3D(T eps2_, T Q_, T invh2_, T dt_, T M_, T h_, T cosTheta_, bool SM_, bool useG_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual solveSBMContactAngle3D<T> * clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    



private:

    T eps2;
    T Q;
    T invh2;
    T dt;
    T M;
    T h;
    T cosTheta;
    bool SM;
    bool useG;

};



#endif // SOLVESBMPROCESSORCONTACTANGLE3D_H_INCLUDED
