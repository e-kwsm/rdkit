#pragma once
#include <memory>
#include <string>

namespace RDKit {
class RWMol;
namespace v2 {
namespace FileParsers {
struct RDKIT_FILEPARSERS_EXPORT CJSONParserParams{};
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromCJSONDataStream(
    std::istream &inStream, unsigned int &line,
    const CJSONParserParams &params = {});
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromCJSONBlock(
    const std::string &molBlock, const CJSONParserParams &params = {});
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromCJSONFile(
    const std::string &fName, const CJSONParserParams &params = {});
}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
