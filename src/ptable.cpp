#include "../include/ptable.h"

// https://en.wikipedia.org/wiki/Atomic_radius#Calculated_atomic_radius
// https://en.wikipedia.org/wiki/Covalent_radius#Radii_for_multiple_bonds
// https://sciencenotes.org/wp-content/uploads/2019/07/CPK-JmolPeriodicTable.pdf
std::unordered_map<std::string, Atom> ptable = {
    {  "H", { .radius = 053.0f, .covalent = 032.0f, .color = { 255.0f / 255.0f, 255.0f / 255.0f, 255.0f / 255.0f }, .mass = 001.0078 } },
    { "Li", { .radius = 167.0f, .covalent = 133.0f, .color = { 204.0f / 255.0f, 128.0f / 255.0f, 255.0f / 255.0f }, .mass = 006.9410 } },
    { "Be", { .radius = 112.0f, .covalent = 102.0f, .color = { 194.0f / 255.0f, 255.0f / 255.0f, 000.0f / 255.0f }, .mass = 009.0122 } },
    {  "C", { .radius = 067.0f, .covalent = 075.0f, .color = { 144.0f / 255.0f, 144.0f / 255.0f, 144.0f / 255.0f }, .mass = 012.0110 } },
    {  "N", { .radius = 056.0f, .covalent = 071.0f, .color = { 048.0f / 255.0f, 080.0f / 255.0f, 248.0f / 255.0f }, .mass = 014.0067 } },
    {  "O", { .radius = 048.0f, .covalent = 063.0f, .color = { 255.0f / 255.0f, 013.0f / 255.0f, 013.0f / 255.0f }, .mass = 015.9994 } },
    {  "P", { .radius = 098.0f, .covalent = 111.0f, .color = { 255.0f / 255.0f, 128.0f / 255.0f, 000.0f / 255.0f }, .mass = 030.9738 } },
    {  "S", { .radius = 088.0f, .covalent = 103.0f, .color = { 255.0f / 255.0f, 255.0f / 255.0f, 048.0f / 255.0f }, .mass = 032.0650 } },
    { "Xe", { .radius = 108.0f, .covalent = 131.0f, .color = { 066.0f / 255.0f, 158.0f / 255.0f, 176.0f / 255.0f }, .mass = 131.2930 } },
    { "El", { .radius = 020.0f, .covalent = 020.0f, .color = { 066.0f / 255.0f, 158.0f / 255.0f, 176.0f / 255.0f }, .mass = 000.0000 } }
};
