#ifndef FILE_FACTORY
#define FILE_FACTORY

#include "Frame.hpp"
#include "ConfigHandler.hpp"
#include "Mesh.hpp"

class FileFactory {
public:
    static void makeFrameFile(Frame* frame, ConfigHandler* config);

    static void makeFile(Frame* frame, ConfigHandler* config, Mesh* mesh, OutputMode mode);
};

#endif /* end of include guard: FILE_FACTORY */
