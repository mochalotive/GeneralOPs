#ifndef H_KERNELPARAMETERS
#define H_KERNELPARAMETERS

#include "crystdist/include/statparameters.h"
#include "converters/include/converters.h"
#include <vector>
#include <string>

enum KernelType {EMPTY_KERNEL, DOPS_KERNEL, BOPS_KERNEL, ROPS_KERNEL};
class KernelParameters {
private:
KernelType kernelType;

public:
CrystalDistributionParameters *distparameters;
KernelParameters(){kernelType=EMPTY_KERNEL;};
KernelParameters(CrystalDistributionParameters *distparameters_, KernelType kernelType_){kernelType = kernelType_, distparameters = distparameters_; };
KernelType returnKernelType(){return(kernelType);};
};

namespace KernelParametersConverter {
//remeber, we're getting passed ALL the strings, and we're appending to an
//array, so keep good track of your input data
//void converter(std::vector<std::string> &s, void *data) {
//std::vector<KernelParameters>* returnKernelParameters = static_cast<std::vector<BasicType>* >(data);
//return;
//  }
//}
//REGISTER_CLASS(KernelParameters, "KERNEL_PARAMETERS",&KernelParametersConverter::converter)
}


#endif
