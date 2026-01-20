#include "plumedinput/include/plumedparameters.h"

std::ostream &operator << ( std::ostream& strm, PlumedFileType tt )
{
   const std::string nameTT[] = { "DOPS", "ROPS", "BOPS", "LOPS", "NOPS", "TOPS"};
   return strm << nameTT[tt];
}

PlumedSystemParameters::PlumedSystemParameters(
                        std::vector<bool> COM_, 
                        std::vector<Real> box_size_, //box size
                        std::vector<std::string> residueNames_,
                        std::vector<Real> cutoff_,
                        std::vector<bool> contact_,
                        std::vector<Real> contact_cutoff_,
                        std::vector<Real> kappa_
                        ):
                        COM(COM_), 
                        box_size(box_size_), 
                        residueNames(residueNames_),
                        cutoff(cutoff_),
                        contact(contact_),
                        contact_cutoff(contact_cutoff_),
                        kappa(kappa_){}




SystemDataPackage::SystemDataPackage(System<Geometry>& mySystem_,
                    std::vector<CrystalDistributionParameters>& myCrystalDistributionParameters_,
                    std::vector<std::vector<MoleculeMap> >& myMoleculeMaps_,
                    std::vector<CrystalOrderParameters>& myCrystalOrderParameters_,
                    std::vector<CrystalOrderParameters>& myOldCrystalOrderParameters_,
                    std::vector<EquationContainer>& myEq_,
                    std::pair<int,int> steps_
                    ):
                    mySystem(mySystem_),
                    myCrystalDistributionParameters(myCrystalDistributionParameters_),
                    myMoleculeMaps(myMoleculeMaps_),
                    myCrystalOrderParameters(myCrystalOrderParameters_),
                    myOldCrystalOrderParameters(myOldCrystalOrderParameters_),
                    myEq(myEq_),
                    steps(steps_)
                   {}

PlumedParameters::PlumedParameters(std::vector<PlumedFileType> fileType_, SystemDataPackage& mySystemDataPackage_, PlumedSystemParameters& myPlumedSystemParameters_, size_t op_, size_t resid_):
                                   fileType(fileType_), mySystemDataPackage(mySystemDataPackage_), myPlumedSystemParameters(myPlumedSystemParameters_), op(op_), resid(resid_)
{
  //system (pointer to System<Geometry>) 
  //activeRes (active residues (compared to all residues), replaced w/ residueNames


  //figures out ALL residue names  
  for (size_t i=0; i< myPlumedSystemParameters.residueNames.size(); i++)
  {
  allRes.push_back(myPlumedSystemParameters.residueNames[i]);
  }
  for (size_t i=0; i<mySystemDataPackage.mySystem.size(); i++)
  {
    bool found=false;
    for (size_t j=0; j<allRes.size(); j++)
    {
      if (compareInternal(allRes[j], mySystemDataPackage.mySystem[i].name()) == 0) found=true;
    }
    if (!found) allRes.push_back(mySystemDataPackage.mySystem[i].name());
  }

  //figures out the numbers of each, regardless of whether we're using them.
  bool molSize;
  for (size_t i=0; i < allRes.size(); i++) 
  {
    molSize=false;
    residueNumbers.push_back(0);
    for(size_t j=0; j<mySystemDataPackage.mySystem.size(); j++)
      {
        if (compareInternal( allRes[i], mySystemDataPackage.mySystem[j].name()  ) == 0 && !molSize) {residueNumbers[i]++; molSize=true;}
        else if(compareInternal( allRes[i], mySystemDataPackage.mySystem[j].name()  ) == 0) residueNumbers[i]++;
      }
  } 
//now know the size of each type of molecule, as well as how many of each


}

int PlumedParameters::compareInternal(std::string s, std::string t)//notice t is a copy :)
  {
   std::string::iterator end_pos = std::remove(t.begin(), t.end(), ' ');
   t.erase(end_pos, t.end()); 
   return (s.compare(t));
  }

PlumedParameters::PlumedParameters(std::vector<std::string> fileType_, SystemDataPackage& mySystemDataPackage_, PlumedSystemParameters& myPlumedSystemParameters_, size_t op_, size_t resid_):
                                   mySystemDataPackage(mySystemDataPackage_), myPlumedSystemParameters(myPlumedSystemParameters_), op(op_), resid(resid_)

{
  
  for (size_t i=0; i<fileType_.size(); i++)
    {
    if (compareInternal(fileType_[i], "DOPS")==0) fileType.push_back(DOPS);
    else if (compareInternal(fileType_[i], "BOPS")==0) fileType.push_back(BOPS);
    else if (compareInternal(fileType_[i], "ROPS")==0) fileType.push_back(ROPS);
    else if (compareInternal(fileType_[i], "LOPS")==0) fileType.push_back(LOPS);
    else if (compareInternal(fileType_[i], "NOPS")==0) fileType.push_back(NOPS);
    else if (compareInternal(fileType_[i], "TOPS")==0) fileType.push_back(TOPS);

    else {std::cout << fileType_[i] << " is not a valid op type. expecting DOPS, ROPS, BOPS, LOPS, NOPS, TOPS." << std::endl; throw std::exception();}
    }

  //figures out ALL residue names  
  for (size_t i=0; i< myPlumedSystemParameters.residueNames.size(); i++)
  {
  allRes.push_back(myPlumedSystemParameters.residueNames[i]);
  }
  for (size_t i=0; i<mySystemDataPackage.mySystem.size(); i++)
  {
    bool found=false;
    for (size_t j=0; j<allRes.size(); j++)
    {
      if (compareInternal(allRes[j], mySystemDataPackage.mySystem[i].name()) == 0) found=true;
    }
    if (!found) allRes.push_back(mySystemDataPackage.mySystem[i].name());
  }

  //figures out the numbers of each, regardless of whether we're using them.
  bool molSize;
  for (size_t i=0; i < allRes.size(); i++) 
  {
    molSize=false;
    residueNumbers.push_back(0);
    for(size_t j=0; j<mySystemDataPackage.mySystem.size(); j++)
      {
        if (compareInternal( allRes[i], mySystemDataPackage.mySystem[j].name()  ) == 0 && !molSize) {residueNumbers[i]++; molSize=true;}
        else if(compareInternal( allRes[i], mySystemDataPackage.mySystem[j].name()  ) == 0) residueNumbers[i]++;
      }
  } 
//now know the size of each type of molecule, as well as how many of each




}
