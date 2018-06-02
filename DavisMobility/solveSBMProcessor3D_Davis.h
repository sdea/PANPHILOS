#ifndef SOLVESBMPROCESSOR3D_DAVIS_DM_H_INCLUDED
#define SOLVESBMPROCESSOR3D_DAVIS_DM_H_INCLUDED

#include "palabos3D.h"
#include "palabos3D.hh"
#include "math.h"

using namespace plb;

template<typename T>
class solveSBM3D_Davis_dM : public BoxProcessingFunctional3D {

public:

    // Constructor
    solveSBM3D_Davis_dM(T eps2_, T Q_, T invh2_, T dt_, T M_, T h_, T cosTheta_, T a1_, T a2_, T a3_, bool SM_, bool useG_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual solveSBM3D_Davis_dM<T> * clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    



private:

    T eps2;
    T Q;
    T invh2;
    T dt;
    T M;
    T h;
    T cosTheta;
    T a1;
    T a2;
    T a3;
    bool SM;
    bool useG;

};



#endif // SOLVESBMPROCESSOR3D_DAVIS_DM_H_INCLUDED
