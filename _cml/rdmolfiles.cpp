#ifndef __clang__version__                     // FIXME
#pragma GCC diagnostic push                    // FIXME
#pragma GCC diagnostic ignored "-Wall"         // FIXME
#pragma GCC diagnostic ignored "-Wextra"       // FIXME
#pragma GCC diagnostic ignored "-Wpedantic"    // FIXME
#else                                          // FIXME
#pragma clang diagnostic push                  // FIXME
#pragma clang diagnostic ignored "-Wall"       // FIXME
#pragma clang diagnostic ignored "-Wextra"     // FIXME
#pragma clang diagnostic ignored "-Wpedantic"  // FIXME
#endif                                         // FIXME

#include <RDBoost/python.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/FileParseException.h>

#include <RDBoost/Wrap.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/BadFileException.h>

#ifndef __clang__version__    // FIXME
#pragma GCC diagnostic pop    // FIXME
#else                         // FIXME
#pragma clang diagnostic pop  // FIXME
#endif                        // FIXME

namespace python = boost::python;
using namespace RDKit;

namespace RDKit {
ROMol* MolFromCMLBlock(const std::string& block, bool sanitize = true,
                       bool removeHs = true) {
  RDLog::InitLogs();  // FIXME
  std::unique_ptr<RWMol> mol;
  v2::FileParsers::CMLFileParserParams params{};
  params.sanitize = sanitize;
  params.removeHs = removeHs;
  try {
    mol = v2::FileParsers::MolFromCMLBlock(block, params);
  } catch (const RDKit::FileParseException& e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol*>(mol.release());
}

ROMol* MolFromCMLFile(const std::string& filename, bool sanitize = true,
                      bool removeHs = true) {
  RDLog::InitLogs();  // FIXME
  std::unique_ptr<RWMol> mol;
  v2::FileParsers::CMLFileParserParams params{};
  params.sanitize = sanitize;
  params.removeHs = removeHs;
  try {
    mol = v2::FileParsers::MolFromCMLFile(filename, params);
  } catch (const RDKit::BadFileException& e) {
    PyErr_SetString(PyExc_IOError, e.what());
    throw boost::python::error_already_set();
  } catch (const RDKit::FileParseException& e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol*>(mol.release());
}
}  // namespace RDKit

BOOST_PYTHON_MODULE(rdmolfiles) {
  std::string docString;

  docString =
      "Construct a molecule from a CML block\n\n\
  ARGUMENTS:\n\
\n\
    - block: string containing CML block\n\
    - sanitize: FIXME\n\
    - removeHs: FIXME\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";
  python::def("MolFromCMLBlock", RDKit::MolFromCMLBlock,
              (python::arg{"block"}, python::arg{"sanitize"} = true,
               python::arg{"removeHs"} = true),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Construct a molecule from a CML file\n\n\
  ARGUMENTS:\n\
\n\
    - filename: name of the file to read\n\
    - sanitize: FIXME\n\
    - removeHs: FIXME\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";
  python::def("MolFromCMLFile", RDKit::MolFromCMLFile,
              (python::arg{"filename"}, python::arg{"sanitize"} = true,
               python::arg{"removeHs"} = true),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
}
