#include "distributor.h"

int main(int argc, char** argv) {
    libint2::initialize(); Distributor(argc, argv).run(); libint2::finalize();
}
