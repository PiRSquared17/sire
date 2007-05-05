// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "SwitchFuncBase.pypp.hpp"
#include "boost/python.hpp"
#include "siremm_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/molecules.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/atom.h"

namespace bp = boost::python;

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireMM::SwitchFuncBase&){ return "SireMM::SwitchFuncBase";}

void register_SwitchFuncBase_class(){

    bp::class_< SireMM::SwitchFuncBase, bp::bases< SireBase::PropertyBase >, boost::noncopyable >( "SwitchFuncBase", bp::no_init )    
        .def( 
            "cutoffDistance"
            , (double ( ::SireMM::SwitchFuncBase::* )(  ) const)( &::SireMM::SwitchFuncBase::cutoffDistance ) )    
        .def( 
            "electrostaticScaleFactor"
            , (double ( ::SireMM::SwitchFuncBase::* )( double ) const)( &::SireMM::SwitchFuncBase::electrostaticScaleFactor )
            , ( bp::arg("dist") ) )    
        .def( 
            "vdwScaleFactor"
            , (double ( ::SireMM::SwitchFuncBase::* )( double ) const)( &::SireMM::SwitchFuncBase::vdwScaleFactor )
            , ( bp::arg("dist") ) )    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMM::SwitchFuncBase >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMM::SwitchFuncBase >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
