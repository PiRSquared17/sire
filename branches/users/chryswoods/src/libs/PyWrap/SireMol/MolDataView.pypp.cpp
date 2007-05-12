// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "MolDataView.pypp.hpp"
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

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireMol::MolDataView&){ return "SireMol::MolDataView";}

void register_MolDataView_class(){

    bp::class_< SireMol::MolDataView, boost::noncopyable >( "MolDataView", bp::no_init )    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMol::MolDataView >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMol::MolDataView >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
