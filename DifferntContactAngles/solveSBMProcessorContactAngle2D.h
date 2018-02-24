#ifndef SOLVESBMPROCESSORCONTACTANGLE2D_H_INCLUDED
#define SOLVESBMPROCESSORCONTACTANGLE2D_H_INCLUDED

#include "palabos2D.h"
#include "palabos2D.hh"

using namespace plb;

template<typename T>
class solveSBMContactAngle2D : public BoxProcessingFunctional2D {

public:

    // Constructor
    solveSBMContactAngle2D(T eps2_, T Q_, T invh2_, T dt_, T M_, T h_, T cosTheta_, bool SM_, bool useG_);
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual solveSBMContactAngle2D<T> * clone() const;
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



#endif // SOLVESBMPROCESSORCONTACTANGLE2D_H_INCLUDED
