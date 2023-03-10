#include "../include/logger.h"

void Logger::Log(bool silent, Eigen::MatrixXd matrix) {
    if (silent) return;
    std::string fmt = "%20.14f";
    std::string format1 = fmt, format2 = matrix.cols() % 5 ? fmt : "";
    for (int i = 0; i < std::min((int)matrix.cols(), 4); i++) format1 += " " + fmt;
    for (int i = 0; i < matrix.cols() % 5 - 1; i++) format2 += " " + fmt;

    for (int i = 0; i <= matrix.cols() / 5; i++) {
        Logger::Log(silent, "%i-%i", i * 5 + 1, std::min(i * 5 + 5, (int)matrix.cols()));
        for (int j = 0; j < matrix.rows(); j++) {
            if (5 * i + 4 < matrix.cols()) Logger::Log(false, format1, matrix(j, 5 * i), matrix(j, 5 * i + 1), matrix(j, 5 * i + 2), matrix(j, 5 * i + 3), matrix(j, 5 * i + 4));
            else if (matrix.cols() % 5 == 1) Logger::Log(false, format2, matrix(j, 5 * i));
            else if (matrix.cols() % 5 == 2) Logger::Log(false, format2, matrix(j, 5 * i), matrix(j, 5 * i + 1));
            else if (matrix.cols() % 5 == 3) Logger::Log(false, format2, matrix(j, 5 * i), matrix(j, 5 * i + 1), matrix(j, 5 * i + 2));
            else if (matrix.cols() % 5 == 4) Logger::Log(false, format2, matrix(j, 5 * i), matrix(j, 5 * i + 1), matrix(j, 5 * i + 2), matrix(j, 5 * i + 3));
        }
    }
}
