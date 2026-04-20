#include <iostream>
#include <cstdlib>
#include "cCubeModel.h"

int main(int argc, char* argv[])
{
    const char* config = (argc > 1) ? argv[1] : "cli/config_cube.xml";

    cCubeModel model;
    model.LoadConfig(config);
    model.Run();

    return EXIT_SUCCESS;
}
