#ifndef SOLVESBMPROCESSOR2D_IMPLICIT_H_INCLUDED
#define SOLVESBMPROCESSOR2D_IMPLICIT_H_INCLUDED

#include "palabos2D.h"
#include "palabos2D.hh"

using namespace plb;

template<typename T>
class solveSBM2DImplicit : public BoxProcessingFunctional2D {

public:

    // Constructor
    solveSBM2DImplicit(T eps2_, T Q_, T invh2_, T dt_, T M_, T h_, T cosTheta_, bool SM_, bool useG_);
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual solveSBM2DImplicit<T> * clone() const;
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



#endif // SOLVESBMPROCESSOR2D_IMPLICIT_H_INCLUDED
