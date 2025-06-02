#include <iosfwd>

namespace RDKit {
class ROMol;

void MolToCJSONBlock(std::ostream &os, const ROMol &mol, int confId,
                     bool kekulize);
std::string MolToCJSONBlock(const ROMol &mol, int confId, bool kekulize);
void MolToCJSONFile(const ROMol &mol, const std::string &fName, int confId,
                    bool kekulize);
}  // namespace RDKit
