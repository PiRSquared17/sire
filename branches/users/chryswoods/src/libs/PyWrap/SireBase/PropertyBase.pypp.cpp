// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "PropertyBase.pypp.hpp"
#include "boost/python.hpp"
#include "sirebase_headers.h"

namespace bp = boost::python;

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireBase::PropertyBase&){ return "SireBase::PropertyBase";}

void register_PropertyBase_class(){

    bp::class_< SireBase::PropertyBase, std::auto_ptr<SireBase::PropertyBase>, boost::noncopyable >( "PropertyBase", bp::no_init )    
        .def( 
            "null_property"
            , (::SireBase::Property (*)(  ))( &::SireBase::PropertyBase::null_property ) )    
        .def( bp::self != bp::self )    
        .def( bp::self == bp::self )    
        .def( 
            "typeName"
            , (char const * (*)(  ))( &::SireBase::PropertyBase::typeName ) )    
        .def( 
            "what"
            , (char const * ( ::SireBase::PropertyBase::* )(  ) const)( &::SireBase::PropertyBase::what ) )    
        .staticmethod( "null_property" )    
        .staticmethod( "typeName" )    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireBase::PropertyBase >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireBase::PropertyBase >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
