// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "DataPoint.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireMaths/maths.h"

#include "SireStream/shareddatastream.h"

#include "fep.h"

#include "tostring.h"

#include "fep.h"

Soiree::DataPoint __copy__(const Soiree::DataPoint &other){ return Soiree::DataPoint(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_DataPoint_class(){

    { //::Soiree::DataPoint
        typedef bp::class_< Soiree::DataPoint > DataPoint_exposer_t;
        DataPoint_exposer_t DataPoint_exposer = DataPoint_exposer_t( "DataPoint", bp::init< >() );
        bp::scope DataPoint_scope( DataPoint_exposer );
        DataPoint_exposer.def( bp::init< double, double >(( bp::arg("x"), bp::arg("y") )) );
        DataPoint_exposer.def( bp::init< double, double, double, double >(( bp::arg("x"), bp::arg("y"), bp::arg("xerror"), bp::arg("yerror") )) );
        DataPoint_exposer.def( bp::init< double, double, double, double, double, double >(( bp::arg("x"), bp::arg("y"), bp::arg("xminerror"), bp::arg("yminerror"), bp::arg("xmaxerror"), bp::arg("ymaxerror") )) );
        DataPoint_exposer.def( bp::init< Soiree::DataPoint const & >(( bp::arg("other") )) );
        { //::Soiree::DataPoint::equalWithinError
        
            typedef bool ( ::Soiree::DataPoint::*equalWithinError_function_type )( ::Soiree::DataPoint const & ) const;
            equalWithinError_function_type equalWithinError_function_value( &::Soiree::DataPoint::equalWithinError );
            
            DataPoint_exposer.def( 
                "equalWithinError"
                , equalWithinError_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::Soiree::DataPoint::equalWithinMaxError
        
            typedef bool ( ::Soiree::DataPoint::*equalWithinMaxError_function_type )( ::Soiree::DataPoint const & ) const;
            equalWithinMaxError_function_type equalWithinMaxError_function_value( &::Soiree::DataPoint::equalWithinMaxError );
            
            DataPoint_exposer.def( 
                "equalWithinMaxError"
                , equalWithinMaxError_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::Soiree::DataPoint::equalWithinMinError
        
            typedef bool ( ::Soiree::DataPoint::*equalWithinMinError_function_type )( ::Soiree::DataPoint const & ) const;
            equalWithinMinError_function_type equalWithinMinError_function_value( &::Soiree::DataPoint::equalWithinMinError );
            
            DataPoint_exposer.def( 
                "equalWithinMinError"
                , equalWithinMinError_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::Soiree::DataPoint::hasError
        
            typedef bool ( ::Soiree::DataPoint::*hasError_function_type )(  ) const;
            hasError_function_type hasError_function_value( &::Soiree::DataPoint::hasError );
            
            DataPoint_exposer.def( 
                "hasError"
                , hasError_function_value );
        
        }
        { //::Soiree::DataPoint::hasErrorRange
        
            typedef bool ( ::Soiree::DataPoint::*hasErrorRange_function_type )(  ) const;
            hasErrorRange_function_type hasErrorRange_function_value( &::Soiree::DataPoint::hasErrorRange );
            
            DataPoint_exposer.def( 
                "hasErrorRange"
                , hasErrorRange_function_value );
        
        }
        { //::Soiree::DataPoint::hasXError
        
            typedef bool ( ::Soiree::DataPoint::*hasXError_function_type )(  ) const;
            hasXError_function_type hasXError_function_value( &::Soiree::DataPoint::hasXError );
            
            DataPoint_exposer.def( 
                "hasXError"
                , hasXError_function_value );
        
        }
        { //::Soiree::DataPoint::hasXErrorRange
        
            typedef bool ( ::Soiree::DataPoint::*hasXErrorRange_function_type )(  ) const;
            hasXErrorRange_function_type hasXErrorRange_function_value( &::Soiree::DataPoint::hasXErrorRange );
            
            DataPoint_exposer.def( 
                "hasXErrorRange"
                , hasXErrorRange_function_value );
        
        }
        { //::Soiree::DataPoint::hasYError
        
            typedef bool ( ::Soiree::DataPoint::*hasYError_function_type )(  ) const;
            hasYError_function_type hasYError_function_value( &::Soiree::DataPoint::hasYError );
            
            DataPoint_exposer.def( 
                "hasYError"
                , hasYError_function_value );
        
        }
        { //::Soiree::DataPoint::hasYErrorRange
        
            typedef bool ( ::Soiree::DataPoint::*hasYErrorRange_function_type )(  ) const;
            hasYErrorRange_function_type hasYErrorRange_function_value( &::Soiree::DataPoint::hasYErrorRange );
            
            DataPoint_exposer.def( 
                "hasYErrorRange"
                , hasYErrorRange_function_value );
        
        }
        DataPoint_exposer.def( bp::self != bp::self );
        { //::Soiree::DataPoint::operator=
        
            typedef ::Soiree::DataPoint & ( ::Soiree::DataPoint::*assign_function_type )( ::Soiree::DataPoint const & ) ;
            assign_function_type assign_function_value( &::Soiree::DataPoint::operator= );
            
            DataPoint_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        DataPoint_exposer.def( bp::self == bp::self );
        { //::Soiree::DataPoint::toString
        
            typedef ::QString ( ::Soiree::DataPoint::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::Soiree::DataPoint::toString );
            
            DataPoint_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::Soiree::DataPoint::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::Soiree::DataPoint::typeName );
            
            DataPoint_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::Soiree::DataPoint::what
        
            typedef char const * ( ::Soiree::DataPoint::*what_function_type )(  ) const;
            what_function_type what_function_value( &::Soiree::DataPoint::what );
            
            DataPoint_exposer.def( 
                "what"
                , what_function_value );
        
        }
        { //::Soiree::DataPoint::x
        
            typedef double ( ::Soiree::DataPoint::*x_function_type )(  ) const;
            x_function_type x_function_value( &::Soiree::DataPoint::x );
            
            DataPoint_exposer.def( 
                "x"
                , x_function_value );
        
        }
        { //::Soiree::DataPoint::xError
        
            typedef double ( ::Soiree::DataPoint::*xError_function_type )(  ) const;
            xError_function_type xError_function_value( &::Soiree::DataPoint::xError );
            
            DataPoint_exposer.def( 
                "xError"
                , xError_function_value );
        
        }
        { //::Soiree::DataPoint::xMaxError
        
            typedef double ( ::Soiree::DataPoint::*xMaxError_function_type )(  ) const;
            xMaxError_function_type xMaxError_function_value( &::Soiree::DataPoint::xMaxError );
            
            DataPoint_exposer.def( 
                "xMaxError"
                , xMaxError_function_value );
        
        }
        { //::Soiree::DataPoint::xMinError
        
            typedef double ( ::Soiree::DataPoint::*xMinError_function_type )(  ) const;
            xMinError_function_type xMinError_function_value( &::Soiree::DataPoint::xMinError );
            
            DataPoint_exposer.def( 
                "xMinError"
                , xMinError_function_value );
        
        }
        { //::Soiree::DataPoint::y
        
            typedef double ( ::Soiree::DataPoint::*y_function_type )(  ) const;
            y_function_type y_function_value( &::Soiree::DataPoint::y );
            
            DataPoint_exposer.def( 
                "y"
                , y_function_value );
        
        }
        { //::Soiree::DataPoint::yError
        
            typedef double ( ::Soiree::DataPoint::*yError_function_type )(  ) const;
            yError_function_type yError_function_value( &::Soiree::DataPoint::yError );
            
            DataPoint_exposer.def( 
                "yError"
                , yError_function_value );
        
        }
        { //::Soiree::DataPoint::yMaxError
        
            typedef double ( ::Soiree::DataPoint::*yMaxError_function_type )(  ) const;
            yMaxError_function_type yMaxError_function_value( &::Soiree::DataPoint::yMaxError );
            
            DataPoint_exposer.def( 
                "yMaxError"
                , yMaxError_function_value );
        
        }
        { //::Soiree::DataPoint::yMinError
        
            typedef double ( ::Soiree::DataPoint::*yMinError_function_type )(  ) const;
            yMinError_function_type yMinError_function_value( &::Soiree::DataPoint::yMinError );
            
            DataPoint_exposer.def( 
                "yMinError"
                , yMinError_function_value );
        
        }
        DataPoint_exposer.staticmethod( "typeName" );
        DataPoint_exposer.def( "__copy__", &__copy__);
        DataPoint_exposer.def( "__deepcopy__", &__copy__);
        DataPoint_exposer.def( "clone", &__copy__);
        DataPoint_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::Soiree::DataPoint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        DataPoint_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::Soiree::DataPoint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        DataPoint_exposer.def( "__str__", &__str__< ::Soiree::DataPoint > );
        DataPoint_exposer.def( "__repr__", &__str__< ::Soiree::DataPoint > );
    }

}