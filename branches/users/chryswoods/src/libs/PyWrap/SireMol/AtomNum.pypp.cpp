// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "AtomNum.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMol/moleculedata.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"
#include "SireVol/space.h"

namespace bp = boost::python;

SireMol::AtomNum __copy__(const SireMol::AtomNum &other){ return SireMol::AtomNum(other); }

#include "SirePy/str.hpp"

void register_AtomNum_class(){

    bp::class_< SireMol::AtomNum, bp::bases< SireMol::IDBase > >( "AtomNum" )    
        .def( bp::init< >() )    
        .def( bp::init< quint32 >(( bp::arg("id") )) )    
        .def( "__copy__", &__copy__)    
        .def( "__str__", &SirePy::__str__< ::SireMol::AtomNum > );

}
