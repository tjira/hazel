#pragma once

#include <boost/algorithm/string.hpp>

#define EH2EV 27.211324570273

#define NELECTRONS(SYS) ([](auto s){auto a=s.getAtoms();return std::accumulate(a.begin(),a.end(),0,[](int e,auto a){return e+a.atomic_number;});}(SYS))
#define NUCREPALL(SYS) ([](auto s){double v=0;for(auto a:s.getAtoms())for(auto b:s.getAtoms())if(a!=b)v+=NUCREP(a, b)/2;return v;}(SYS))
#define NUCREP(A, B) (a.atomic_number*b.atomic_number/Eigen::Vector3d(a.x-b.x,a.y-b.y,a.z-b.z).norm())

#define FCONTENTS(F) ([](auto f){std::ifstream s(f);std::stringstream b;b<<s.rdbuf();auto o=b.str();boost::trim_right(o);return o;}(F))
#define PABSOLUTE(F) (std::filesystem::absolute(argv[1]).parent_path().string()+"/"+F)

#define STRINGIFY(X) STRING(X)
#define STRING(X) #X
