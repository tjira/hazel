#include "../include/libint.h"

bool libint2::operator!=(const libint2::Atom& a, const libint2::Atom& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.atomic_number != b.atomic_number;
}
