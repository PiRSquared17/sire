// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "VariantProperty.pypp.hpp"
#include "boost/python.hpp"
#include "sirebase_headers.h"

namespace bp = boost::python;

SireBase::VariantProperty __copy__(const SireBase::VariantProperty &other){ return SireBase::VariantProperty(other); }

#include "SireQt/qdatastream.hpp"

#include "SirePy/str.hpp"

void register_VariantProperty_class(){

    bp::class_< SireBase::VariantProperty, bp::bases< QVariant, SireBase::ConcreteProperty<SireBase::VariantProperty, SireBase::PropertyBase>, SireBase::PropertyBase > >( "VariantProperty" )    
        .def( bp::init< >() )    
        .def( bp::init< QVariant const & >(( bp::arg("value") )) )    
        .def( bp::init< SireBase::Property const & >(( bp::arg("other") )) )    
        .def( bp::self != bp::self )    
        .def( bp::self == bp::self )    
        .def( 
            "typeName"
            , (char const * (*)(  ))( &::SireBase::VariantProperty::typeName ) )    
        .staticmethod( "typeName" )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireBase::VariantProperty >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireBase::VariantProperty >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &SirePy::__str__< ::SireBase::VariantProperty > );

}
