#include "plumedinput/include/kernelfile.h"
KernelFile::KernelFile():IOFile(){};
KernelFile::KernelFile(std::string const &fileName):IOFile(fileName, OUT){};

KernelFile &operator<<(KernelFile &kernelFile, KernelParameters &parameters)
{
  assert(kernelFile.is_open() && kernelFile.mode() == OUT);
  kernelFile.setf(std::ios::right, std::ios::adjustfield);
  std::streamsize precision = kernelFile.precision();
  kernelFile.precision(WRITE_PRECISION);
switch (parameters.returnKernelType()){
  case DOPS_KERNEL:  
    kernelFile << "#! FIELDS height mu sigma"  <<std::endl;
    kernelFile << "#! SET kerneltype gaussian" << std::endl;
    break;
  case BOPS_KERNEL:
    kernelFile << "#! FIELDS height kappa mu_i mu_j mu_k"  <<std::endl;
    kernelFile << "#! SET kerneltype gaussian" << std::endl;
    break;
  case ROPS_KERNEL:
    kernelFile << "#! FIELDS height kappa mu_w mu_i mu_j mu_k"  <<std::endl;
    kernelFile << "#! SET kerneltype gaussian" << std::endl;
    break;
  default:
    std::cerr << "Kernel File has EMPTY_KERNEL flag set. This usually means you've defined an empty kernelfile by accident." << std::endl;
    return(kernelFile);
}
//same switch again, now to write relevant concentration parameters or variances.
size_t numGroups = parameters.distparameters->means.size();//can do this ahead of time 
switch (parameters.returnKernelType()){
  case DOPS_KERNEL:
    {
    for (size_t i=0; i<numGroups; i++){
      size_t numPeaks = parameters.distparameters->means[i].size();
      for (size_t j=0; j<numPeaks; j++){
        kernelFile << parameters.distparameters->normalizationFactors[i][j].distance << " " << parameters.distparameters->means[i][j].distance << " " << 1/sqrt(parameters.distparameters->concentrationParameters[i][j].distance) << std::endl;
      }
    }
    break;
    }
  case BOPS_KERNEL:
    {
    for (size_t i=0; i<numGroups; i++){
      size_t numPeaks = parameters.distparameters->means[i].size();
      for (size_t j=0; j<numPeaks; j++){
        kernelFile << parameters.distparameters->normalizationFactors[i][j].bondOrientation << " " << parameters.distparameters->concentrationParameters[i][j].bondOrientation << " " << parameters.distparameters->means[i][j].bondOrientation.x << " "  << parameters.distparameters->means[i][j].bondOrientation.y << " " << parameters.distparameters->means[i][j].bondOrientation.z   << std::endl;
        }
      }
    break;
    }
  case ROPS_KERNEL:
    {
    for (size_t i=0; i<numGroups; i++){
      size_t numPeaks = parameters.distparameters->means[i].size();
      for (size_t j=0; j<numPeaks; j++){
        kernelFile << parameters.distparameters->normalizationFactors[i][j].relativeOrientation << " " << parameters.distparameters->concentrationParameters[i][j].relativeOrientation << " " << parameters.distparameters->means[i][j].relativeOrientation.w<< " " <<parameters.distparameters->means[i][j].relativeOrientation.x << " "  << parameters.distparameters->means[i][j].relativeOrientation.y << " " << parameters.distparameters->means[i][j].relativeOrientation.z   << std::endl;
    }
      }    
    break;
    }
  
  default:
    std::cerr << "I'm not even sure how you got here. Kernel File has EMPTY_KERNEL flag set. This usually means you've defined an empty kernelfile by accident." << std::endl;
    return(kernelFile);
}



//



return (kernelFile);
}
