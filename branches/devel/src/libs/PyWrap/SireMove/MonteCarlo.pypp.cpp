// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "MonteCarlo.pypp.hpp"
#include "boost/python.hpp"
#include "siremove_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/moleculegroup.h"
#include "SireSystem/simsystem.h"

namespace bp = boost::python;

void register_MonteCarlo_class(){

    { //::SireMove::MonteCarlo
        typedef bp::class_< SireMove::MonteCarlo, bp::bases< SireSystem::MoveBase >, boost::noncopyable > MonteCarlo_exposer_t;
        MonteCarlo_exposer_t MonteCarlo_exposer = MonteCarlo_exposer_t( "MonteCarlo", bp::no_init );
        bp::scope MonteCarlo_scope( MonteCarlo_exposer );
        bp::class_< SireMove::MonteCarlo::CheckPoint >( "CheckPoint" );
        { //::SireMove::MonteCarlo::generator
        
            typedef ::SireMaths::RanGenerator const & ( ::SireMove::MonteCarlo::*generator_function_type )(  ) const;
            
            MonteCarlo_exposer.def( 
                "generator"
                , generator_function_type( &::SireMove::MonteCarlo::generator )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMove::MonteCarlo::setGenerator
        
            typedef void ( ::SireMove::MonteCarlo::*setGenerator_function_type )( ::SireMaths::RanGenerator const & ) ;
            
            MonteCarlo_exposer.def( 
                "setGenerator"
                , setGenerator_function_type( &::SireMove::MonteCarlo::setGenerator )
                , ( bp::arg("generator") ) );
        
        }
    }

}
