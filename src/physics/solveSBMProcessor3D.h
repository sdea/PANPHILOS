#ifndef SOLVESBMPROCESSOR3D_H_INCLUDED
#define SOLVESBMPROCESSOR3D_H_INCLUDED

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace plb;

template<typename T>
class solveSBM3D : public BoxProcessingFunctional3D {

public:

    // Constructor
    solveSBM3D(T eps2_, T Q_, T invh2_, T dt_, T M_, T h_, T cosTheta_, bool SM_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual solveSBM3D<T> * clone() const;
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

};



#endif // SOLVESBMPROCESSOR3D_H_INCLUDED
