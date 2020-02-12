#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "hdf5_lib.hpp"

int main(int argc, char** argv) {
    h5_file file("/home/abiakarn/git/usc-collab-upstream/arbor/network/h5_files/data_network_s1000.h5");
    file.print();
}