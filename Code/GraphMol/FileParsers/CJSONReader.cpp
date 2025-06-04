#include "cjson.h"
#include <GraphMol/RWMol.h>

namespace RDKit {
namespace v2 {
namespace FileParsers {
std::unique_ptr<RWMol> MolFromCJSONDataStream(std::istream &inStream,
                                              unsigned int &line,
                                              const CJSONParserParams &params) {
  return nullptr;
}

std::unique_ptr<RWMol> MolFromCJSONBlock(const std::string &molBlock,
                                         const CJSONParserParams &params) {
  auto rwmol = std::make_unique<RWMol>();
  return rwmol;
}

std::unique_ptr<RWMol> MolFromCJSONFile(const std::string &fName,
                                        const CJSONParserParams &params) {
  return nullptr;
}
}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
