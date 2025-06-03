#include <iosfwd>
#include <RDGeneral/export.h>

namespace RDKit {
class ROMol;
enum class CJSONCoords {
  _3d,            //!< `[x, y, z, x, y, z, …]`
  _3dFractional,  //!< `[x, y, z, x, y, z, …]`, fractional
  _3dSets         //!< `[[x, y, z], [x, y, z], …]`
};

struct RDKIT_FILEPARSERS_EXPORT CJSONWriterParams {
  bool kekulize = true;
  CJSONCoords coords = CJSONCoords::_3dSets;
};

void MolToCJSONBlock(std::ostream &os, const ROMol &mol,
                     const CJSONWriterParams &params, int confId = -1);
std::string MolToCJSONBlock(const ROMol &mol, const CJSONWriterParams &params,
                            int confId = -1);
void MolToCJSONFile(const ROMol &mol, const std::string &fName,
                    const CJSONWriterParams &params, int confId = -1);
}  // namespace RDKit
