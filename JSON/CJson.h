#include <iosfwd>
#include <RDGeneral/export.h>

namespace RDKit {
class ROMol;
struct RDKIT_FILEPARSERS_EXPORT CJSONWriterParams {
  bool kekulize = true;
};

void MolToCJSONBlock(std::ostream &os, const ROMol &mol,
                     const CJSONWriterParams &params, int confId = -1);
std::string MolToCJSONBlock(const ROMol &mol, const CJSONWriterParams &params,
                            int confId = -1);
void MolToCJSONFile(const ROMol &mol, const std::string &fName,
                    const CJSONWriterParams &params, int confId = -1);
}  // namespace RDKit
