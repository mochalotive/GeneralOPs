#ifndef H_PDBFILE
#define H_PDBFILE

#include "common/include/iofile.h"
#include "mymol/include/geometry.h"
#include "mymol/include/molecule.h"
#include "mymol/include/system.h"

class PSFFile: public IOFile
{
public:

 Constructors

   PSFFile();                            // Defines an empty PSFFile
     PSFFile(std::string const &fileName); // Defines a PSF File with a given file name

     // Interface

     void setFile(std::string const &fileName); // Sets the file name

       // Read to Geometry objects

     friend PSFFile &operator>>(PSFFile &psfFile, Geometry &geometry); // Read geometry from file

         // Read to System<Geometry> objects

     friend PSFFile &operator>>(PSFFile &psfFile, System<Geometry> &system); // Read system from file

           // Read to System<Molecule> objects (inheritance not working with template)

     friend PSFFile &operator>>(PSFFile &psfFile, System<Molecule> &system); // Read system from file

     private:

     void setFile(std::string const &fileName, // Format is read-only - blocking
          IOMode const &mode,          // the general form from IOFile
          IOFormat const &format);
                          };

            #endif

          /*
          ** End of class PSFFile
          */

     // Inline//

     inline void PSFFile::setFile(std::string const &fileName)
        {
          IOFile::setFile(fileName, IN);
        }

     inline void PSFFile::setFile(std::string const &fileName,
                                  IOMode const &mode,
                                  IOFormat const &format)
     {
       std::cerr << "Error in PSFFile: Writing not implemented!" << std::endl;
      }

                                                                                                                                             
