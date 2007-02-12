// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "CuttingFunctionBase.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"

namespace bp = boost::python;

void register_CuttingFunctionBase_class(){

    bp::class_< SireMol::CuttingFunctionBase, boost::noncopyable >( "CuttingFunctionBase", bp::no_init )    
        .def( 
            "__call__"
            , (::SireMol::CutGroupNum ( ::SireMol::CuttingFunctionBase::* )( ::SireMol::AtomIndex const &,::SireMol::EditMolData const & ) const)( &::SireMol::CuttingFunctionBase::operator() )
            , ( bp::arg("atom"), bp::arg("moldata") ) )    
        .def( 
            "what"
            , (char const * ( ::SireMol::CuttingFunctionBase::* )(  ) const)( &::SireMol::CuttingFunctionBase::what ) );

}
