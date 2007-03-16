// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "using_parameters_base.pypp.hpp"
#include "boost/python.hpp"
#include "siredb_headers.h"
#include "SireCAS/function.h"
#include "SireCAS/expression.h"

namespace bp = boost::python;

SireDB::using_parameters_base __copy__(const SireDB::using_parameters_base &other){ return SireDB::using_parameters_base(other); }

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireDB::using_parameters_base&){ return "SireDB::using_parameters_base";}

void register_using_parameters_base_class(){

    bp::class_< SireDB::using_parameters_base, bp::bases< SireDB::using_database > >( "using_parameters_base" )    
        .def( bp::init< >() )    
        .def( bp::init< QStringList const & >(( bp::arg("typenames") )) )    
        .def( 
            "what"
            , (char const * ( ::SireDB::using_parameters_base::* )(  ) const)( &::SireDB::using_parameters_base::what ) )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireDB::using_parameters_base >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireDB::using_parameters_base >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
