// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "ConstantPower.pypp.hpp"
#include "boost/python.hpp"
#include "sirecas_headers.h"

namespace bp = boost::python;

void register_ConstantPower_class(){

    bp::class_< SireCAS::ConstantPower, bp::bases< SireCAS::PowerFunction >, boost::noncopyable >( "ConstantPower", bp::no_init )    
        .def( 
            "core"
            , (::SireCAS::Expression ( ::SireCAS::ConstantPower::* )(  ) const)( &::SireCAS::ConstantPower::core ) )    
        .def( 
            "hash"
            , (::uint ( ::SireCAS::ConstantPower::* )(  ) const)( &::SireCAS::ConstantPower::hash ) );

}
