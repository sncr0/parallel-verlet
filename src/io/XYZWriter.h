#ifndef XYZ_WRITER_H
#define XYZ_WRITER_H

#include <fstream>
#include <string>
#include "../system/System.h"

class XYZWriter {
public:
    explicit XYZWriter(const std::string& filename);
    ~XYZWriter();

    void writeHeader(size_t numParticles);
    void writeFrame(const System& system);
    void close();

private:
    std::ofstream file;
    size_t numParticles;
    bool headerWritten;
};

#endif // XYZ_WRITER_H
