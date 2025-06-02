#include <iosfwd>

namespace RDKit {
class ROMol;

void MolToCJSONBlock(std::ostream &os, const ROMol &mol, int confId = -1,
                     bool kekulize = true);
std::string MolToCJSONBlock(const ROMol &mol, int confId = -1,
                            bool kekulize = true);
void MolToCJSONFile(const ROMol &mol, const std::string &fName, int confId = -1,
                    bool kekulize = true);
}  // namespace RDKit
