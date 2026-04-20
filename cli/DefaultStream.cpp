#include <iostream>
#include <fstream>

// Redirect std::clog to a file so that verbose library output does not
// clutter the terminal.  The ofstream is kept open for the lifetime of
// the process via a static local.
struct DefaultStreamInit {
    DefaultStreamInit() {
        static std::ofstream log_file("atom_cube.log", std::ofstream::out);
        std::clog.rdbuf(log_file.rdbuf());
    }
} default_stream_init;
