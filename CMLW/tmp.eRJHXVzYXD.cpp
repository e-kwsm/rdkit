#include "CMLWriter.h"
#include <fstream>

#include <RDGeneral/BadFileException.h>
#include <RDGeneral/FileParseException.h>

namespace RDKit {
CMLWriter::CMLWriter(const std::string &fileName)
    : CMLWriter{std::make_unique<std::ofstream>(fileName)} {}

CMLWriter::CMLWriter(std::unique_ptr<std::ostream> &&p_os)
    : p_ostream{std::move(p_os)} {
  PRECONDITION(p_os, "null stream");
  if (p_os->bad()) {
    throw FileParseException("Bad output stream");
  }
}

CMLWriter::~CMLWriter() { close(); }

void CMLWriter::write(const ROMol &mol, int confId) {
  PRECONDITION(p_ostream, "no output stream");
}

void CMLWriter::flush() const {
  if (p_ostream) {
    p_ostream->flush();
  }
}

void CMLWriter::close() {
  flush();
  p_ostream.reset();
}
}  // namespace RDKit
