#ifndef COMPUTESEGMENTATION2D_H_INCLUDED
#define COMPUTESEGMENTATION2D_H_INCLUDED

namespace plb {

template<typename T>
class ComputeSegmentation2D : public BoxProcessingFunctional2D {

public:

    ComputeSegmentation2D();
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual ComputeSegmentation2D<T> * clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;

private:


};

}
#endif // COMPUTESEGMENTATION2D_H_INCLUDED
