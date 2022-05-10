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

  /// name of the file where the data is read
  const std::string _file_name;

  /// whether the file contains a header with the column names
  const MooseUtils::DelimitedFileReader::HeaderFlag _header;

  /// string used as a delimiter
  const std::string _delimiter;

  std::vector<std::vector<double>> _data;
};

#endif
