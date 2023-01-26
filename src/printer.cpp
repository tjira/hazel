#include "../include/printer.h"

int Printer::w1 = 120, Printer::w2 = 100, Printer::w3 = 80;

void Printer::printTitle() {
    std::cout << "╔" << repeat("═", w1 - 2) << "╗\n";
    std::cout << "║" << std::string((w1 - 22) / 2 - 1, ' ') << "QUANTUM HAZEL SOFTWARE" << std::string((w1 - 22) / 2 - 1, ' ') << "║\n";
    std::cout << "╚" << repeat("═", w1 - 2) << "╝" << std::endl;
}

std::string Printer::repeat(const std::string& str, int i) {
    std::string output;
    for (int j = 0; j < i; j++) {
        output += str;
    }
    return output;
} 
