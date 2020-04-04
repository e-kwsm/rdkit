#pragma once

#include <RDGeneral/export.h>
#include <RDGeneral/types.h>
#include <GraphMol/RDKitBase.h>

namespace RDKit {
RDKIT_FILEPARSERS_EXPORT RDKit::RWMol* CMLBlockToMol(
    const std::string& block, bool sanitize = true,
    bool removeHs = true) noexcept(false);

RDKIT_FILEPARSERS_EXPORT RDKit::RWMol* CMLFileToMol(
    const std::string& filename, bool sanitize = true,
    bool removeHs = true) noexcept(false);
}  // namespace RDKit
