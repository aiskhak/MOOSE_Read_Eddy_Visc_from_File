#include "AuxVarFromCSVFile.h"

registerMooseObject("MooseApp", AuxVarFromCSVFile);

InputParameters
AuxVarFromCSVFile::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredParam<FileName>("file_name", "where eddy viscosity is stored, x, y, z, mu_t");
  params.addParam<bool>("header", "True if there is a header in the file");
  params.addParam<std::string>("delim", ",", "delimiter between the numbers");
  params.ignoreParameter<ExecFlagEnum>("execute_on");      // read only once
  params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL}; // read only once
  return params;
}

AuxVarFromCSVFile::AuxVarFromCSVFile(const InputParameters & parameters)
  : AuxKernel(parameters),
  _file_name(getParam<FileName>("file_name")),
  _header(isParamValid("header")
              ? (getParam<bool>("header") ? MooseUtils::DelimitedFileReader::HeaderFlag::ON
                                          : MooseUtils::DelimitedFileReader::HeaderFlag::OFF)
              : MooseUtils::DelimitedFileReader::HeaderFlag::AUTO),
  _delim(getParam<std::string>("delim"))
{
  MooseUtils::DelimitedFileReader file(_file_name);
  file.setHeaderFlag(_header);
  file.setDelimiter(_delim);
  file.read();
  _data = file.getData();
}

Real
AuxVarFromCSVFile::computeValue()
{
  // find nearest point
  Real x_qp = _q_point[_qp](0);
  Real y_qp = _q_point[_qp](1);
  Real z_qp = _q_point[_qp](2);
  Real min = 1.e12;
  Real val = 0.0;
  for (unsigned int i = 0; i < _data[0].size(); i++)
  {
    Real dx   = x_qp - _data[0][i];
    Real dy   = y_qp - _data[1][i];
    Real dz   = z_qp - _data[2][i];
    Real dist = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (dist < min)
    {
      min = dist;
      val = _data[3][i];
    }
  }

  return val;
}