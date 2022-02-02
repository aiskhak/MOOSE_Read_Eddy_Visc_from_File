#ifndef AUXVARFROMCSVFILE_H
#define AUXVARFROMCSVFILE_H

#include "AuxKernel.h"
#include "DelimitedFileReader.h"

class AuxVarFromCSVFile : public AuxKernel
{
public:
  static InputParameters validParams();

  AuxVarFromCSVFile(const InputParameters & parameters);

protected:
  virtual Real computeValue();
  const std::string _file_name;
  const MooseUtils::DelimitedFileReader::HeaderFlag _header;
  const std::string _delim;
  std::vector<std::vector<double>> _data;
};

#endif
