#include "../include/ptable.h"

// https://en.wikipedia.org/wiki/Atomic_radius#Calculated_atomic_radius
// https://en.wikipedia.org/wiki/Covalent_radius#Radii_for_multiple_bonds
// https://sciencenotes.org/wp-content/uploads/2019/07/CPK-JmolPeriodicTable.pdf
std::unordered_map<std::string, Atom> ptable = {
    { "H", { .radius = 53.0f, .covalent = 032.0f, .color = { 255.0f / 255.0f, 255.0f / 255.0f, 255.0f / 255.0f }, .mass = 01.0078 } },
    { "C", { .radius = 67.0f, .covalent = 075.0f, .color = { 144.0f / 255.0f, 144.0f / 255.0f, 144.0f / 255.0f }, .mass = 12.0110 } },
    { "N", { .radius = 56.0f, .covalent = 071.0f, .color = { 048.0f / 255.0f, 080.0f / 255.0f, 248.0f / 255.0f }, .mass = 14.0067 } },
    { "O", { .radius = 48.0f, .covalent = 063.0f, .color = { 255.0f / 255.0f, 013.0f / 255.0f, 013.0f / 255.0f }, .mass = 15.9994 } },
    { "P", { .radius = 98.0f, .covalent = 111.0f, .color = { 255.0f / 255.0f, 128.0f / 255.0f, 000.0f / 255.0f }, .mass = 30.9738 } },
    { "S", { .radius = 88.0f, .covalent = 103.0f, .color = { 255.0f / 255.0f, 255.0f / 255.0f, 048.0f / 255.0f }, .mass = 32.0650 } }
};
