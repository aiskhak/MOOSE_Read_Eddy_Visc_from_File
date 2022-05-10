#include "AuxVarFromCSVFile.h"

registerMooseObject("MooseApp", AuxVarFromCSVFile);

InputParameters
AuxVarFromCSVFile::validParams()
{
  InputParameters params = AuxKernel::validParams();

  params.addRequiredParam<FileName>("file_name",
                                    "name of the file in which the time sequence is read");
  params.addParam<bool>("header",
                        "indicates whether the file contains a header with the column names");
  params.addParam<std::string>("delimiter", ",", "delimiter used to parse the file");

  // Force this AuxKernel to perform computeValue() on INITIAL, by ignoring any user input
  params.ignoreParameter<ExecFlagEnum>("execute_on");
  params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL};

  return params;
}

AuxVarFromCSVFile::AuxVarFromCSVFile(const InputParameters & parameters)
  : AuxKernel(parameters),
  _file_name(getParam<FileName>("file_name")),
  _header(isParamValid("header")
              ? (getParam<bool>("header") ? MooseUtils::DelimitedFileReader::HeaderFlag::ON
                                          : MooseUtils::DelimitedFileReader::HeaderFlag::OFF)
              : MooseUtils::DelimitedFileReader::HeaderFlag::AUTO),
  _delimiter(getParam<std::string>("delimiter"))
{
  MooseUtils::DelimitedFileReader file(_file_name);

  file.setHeaderFlag(_header);
  file.setDelimiter(_delimiter);
  file.read();

  _data = file.getData();

  if (_data.size() != 4)
    mooseError("AuxVarFromCSVFile expects that data are given as: 'x, y, z, val'.");
}

Real
AuxVarFromCSVFile::computeValue()
{
  // Find the nearest point from _data, and use its value
  Real x_qp = _q_point[_qp](0);
  Real y_qp = _q_point[_qp](1);
  Real z_qp = _q_point[_qp](2);
  Real min_distance = 1.0e30;
  Real val = 0.0;

  for (unsigned int i = 0; i < _data[0].size(); i++)
  {
    Real dx = x_qp - _data[0][i];
    Real dy = y_qp - _data[1][i];
    Real dz = z_qp - _data[2][i];

    Real distance = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (distance < min_distance)
    {
      min_distance = distance;
      val = _data[3][i];
    }
  }

  return val;
}