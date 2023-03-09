#pragma once

#define EH2EV 27.211324570273

#define NELECTRONS(SYS) ([](auto s){auto a=s.getAtoms();return std::accumulate(a.begin(),a.end(),0,[](int e,auto a){return e+a.atomic_number;});}(SYS))
#define NUCREP(SYS) ([](auto s){double v=0;for(auto a:s.getAtoms())for(auto b:s.getAtoms())if(a!=b)v+=NUCREPOP(a, b)/2;return v;}(SYS))
#define NUCREPOP(A, B) (a.atomic_number*b.atomic_number/Eigen::Vector3d(a.x-b.x,a.y-b.y,a.z-b.z).norm())

#define STRINGIFY(X) STRING(X)
#define STRING(X) #X
