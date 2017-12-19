#ifndef COMPUTEMEANCURVATURE_H_INCLUDED
#define COMPUTEMEANCURVATURE_H_INCLUDED

namespace plb {

template<typename T>
class ComputeMeanCurvature : public BoxProcessingFunctional3D {

public:

    ComputeMeanCurvature(T h_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual ComputeMeanCurvature<T> * clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;

private:

    T h;
};

}
#endif // COMPUTEMEANCURVATURE_H_INCLUDED
