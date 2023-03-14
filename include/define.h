#pragma once

#include <boost/algorithm/string.hpp>
#include <unsupported/Eigen/CXX11/Tensor>
#include <filesystem>

#define DATADIR std::filesystem::current_path()

#define EH2EV 27.211324570273
#define WIDTH 104

#define FCONTENTS(F) ([](auto f){std::ifstream s(f);std::stringstream b;b<<s.rdbuf();auto o=b.str();boost::trim_right(o);return o;}(F))
#define PABSOLUTE(F) (std::filesystem::absolute(program.get<std::string>("input")).parent_path().string()+"/"+F)

#define STRINGIFY(X) STRING(X)
#define STRING(X) #X

typedef Eigen::Vector3d Vec3;
typedef Eigen::VectorXd Vec;
typedef Eigen::MatrixXd Mat;
typedef Eigen::Tensor<double, 4> Ten4;
