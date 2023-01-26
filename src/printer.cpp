#include "../include/printer.h"

int Printer::w1 = 120, Printer::w2 = 100, Printer::w3 = 80;

void Printer::printTitle() {
    std::cout << "╔══════════════════════════════════════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                                      QUANTUM HAZEL SOFTWARE                                      ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════════════════════════════════════╝\n";
    std::cout << std::flush;
}

std::string Printer::repeat(const std::string& str, int i) {
    std::string output;
    for (int j = 0; j < i; j++) {
        output += str;
    }
    return output;
} 
