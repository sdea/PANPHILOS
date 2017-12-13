#ifndef COMPUTESEGMENTATION_H_INCLUDED
#define COMPUTESEGMENTATION_H_INCLUDED

namespace plb {

template<typename T>
class ComputeSegmentation : public BoxProcessingFunctional3D {

public:

    ComputeSegmentation();
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual ComputeSegmentation<T> * clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;

private:


};

}
#endif // COMPUTESEGMENTATION_H_INCLUDED
