#pragma once

#include "define.h"
#include <boost/format.hpp>
#include <Eigen/Eigen>
#include <iostream>

struct Logger {
    template <class... Args> static void Log(bool silent, std::string fmt, Args... args);
    static void Log(bool silent, Eigen::MatrixXd matrix, std::string fmt = "%20.14f");
};

template <class... Args>
void Logger::Log(bool silent, std::string fmt, Args... args) {
    if (!silent) std::cout << (boost::format(fmt) % ... % args) << std::endl;
}
