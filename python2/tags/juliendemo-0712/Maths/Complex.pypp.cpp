// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Complex.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "complex.h"

#include "rational.h"

#include <gsl/gsl_complex_math.h>

#include "complex.h"

SireMaths::Complex __copy__(const SireMaths::Complex &other){ return SireMaths::Complex(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_Complex_class(){

    { //::SireMaths::Complex
        typedef bp::class_< SireMaths::Complex > Complex_exposer_t;
        Complex_exposer_t Complex_exposer = Complex_exposer_t( "Complex", bp::init< bp::optional< double, double > >(( bp::arg("r")=0.0, bp::arg("i")=0.0 )) );
        bp::scope Complex_scope( Complex_exposer );
        Complex_exposer.def( bp::init< gsl_complex const & >(( bp::arg("complex") )) );
        Complex_exposer.def( bp::init< SireMaths::Complex const & >(( bp::arg("other") )) );
        { //::SireMaths::Complex::abs
        
            typedef double ( ::SireMaths::Complex::*abs_function_type )(  ) const;
            abs_function_type abs_function_value( &::SireMaths::Complex::abs );
            
            Complex_exposer.def( 
                "abs"
                , abs_function_value );
        
        }
        { //::SireMaths::Complex::abs2
        
            typedef double ( ::SireMaths::Complex::*abs2_function_type )(  ) const;
            abs2_function_type abs2_function_value( &::SireMaths::Complex::abs2 );
            
            Complex_exposer.def( 
                "abs2"
                , abs2_function_value );
        
        }
        { //::SireMaths::Complex::arg
        
            typedef double ( ::SireMaths::Complex::*arg_function_type )(  ) const;
            arg_function_type arg_function_value( &::SireMaths::Complex::arg );
            
            Complex_exposer.def( 
                "arg"
                , arg_function_value );
        
        }
        { //::SireMaths::Complex::conjugate
        
            typedef ::SireMaths::Complex ( ::SireMaths::Complex::*conjugate_function_type )(  ) const;
            conjugate_function_type conjugate_function_value( &::SireMaths::Complex::conjugate );
            
            Complex_exposer.def( 
                "conjugate"
                , conjugate_function_value );
        
        }
        { //::SireMaths::Complex::imag
        
            typedef double ( ::SireMaths::Complex::*imag_function_type )(  ) const;
            imag_function_type imag_function_value( &::SireMaths::Complex::imag );
            
            Complex_exposer.def( 
                "imag"
                , imag_function_value );
        
        }
        { //::SireMaths::Complex::inverse
        
            typedef ::SireMaths::Complex ( ::SireMaths::Complex::*inverse_function_type )(  ) const;
            inverse_function_type inverse_function_value( &::SireMaths::Complex::inverse );
            
            Complex_exposer.def( 
                "inverse"
                , inverse_function_value );
        
        }
        { //::SireMaths::Complex::isPurelyComplex
        
            typedef bool ( ::SireMaths::Complex::*isPurelyComplex_function_type )(  ) const;
            isPurelyComplex_function_type isPurelyComplex_function_value( &::SireMaths::Complex::isPurelyComplex );
            
            Complex_exposer.def( 
                "isPurelyComplex"
                , isPurelyComplex_function_value );
        
        }
        { //::SireMaths::Complex::isReal
        
            typedef bool ( ::SireMaths::Complex::*isReal_function_type )(  ) const;
            isReal_function_type isReal_function_value( &::SireMaths::Complex::isReal );
            
            Complex_exposer.def( 
                "isReal"
                , isReal_function_value );
        
        }
        { //::SireMaths::Complex::isZero
        
            typedef bool ( ::SireMaths::Complex::*isZero_function_type )(  ) const;
            isZero_function_type isZero_function_value( &::SireMaths::Complex::isZero );
            
            Complex_exposer.def( 
                "isZero"
                , isZero_function_value );
        
        }
        { //::SireMaths::Complex::logAbs
        
            typedef double ( ::SireMaths::Complex::*logAbs_function_type )(  ) const;
            logAbs_function_type logAbs_function_value( &::SireMaths::Complex::logAbs );
            
            Complex_exposer.def( 
                "logAbs"
                , logAbs_function_value );
        
        }
        { //::SireMaths::Complex::negative
        
            typedef ::SireMaths::Complex ( ::SireMaths::Complex::*negative_function_type )(  ) const;
            negative_function_type negative_function_value( &::SireMaths::Complex::negative );
            
            Complex_exposer.def( 
                "negative"
                , negative_function_value );
        
        }
        Complex_exposer.def( bp::self != bp::self );
        Complex_exposer.def( bp::self != bp::other< double >() );
        Complex_exposer.def( -bp::self );
        { //::SireMaths::Complex::operator=
        
            typedef ::SireMaths::Complex & ( ::SireMaths::Complex::*assign_function_type )( ::SireMaths::Complex const & ) ;
            assign_function_type assign_function_value( &::SireMaths::Complex::operator= );
            
            Complex_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        { //::SireMaths::Complex::operator=
        
            typedef ::SireMaths::Complex & ( ::SireMaths::Complex::*assign_function_type )( double ) ;
            assign_function_type assign_function_value( &::SireMaths::Complex::operator= );
            
            Complex_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("r") )
                , bp::return_self< >() );
        
        }
        Complex_exposer.def( bp::self == bp::self );
        Complex_exposer.def( bp::self == bp::other< double >() );
        { //::SireMaths::Complex::polar
        
            typedef ::SireMaths::Complex ( *polar_function_type )( double,double );
            polar_function_type polar_function_value( &::SireMaths::Complex::polar );
            
            Complex_exposer.def( 
                "polar"
                , polar_function_value
                , ( bp::arg("r"), bp::arg("theta") ) );
        
        }
        { //::SireMaths::Complex::real
        
            typedef double ( ::SireMaths::Complex::*real_function_type )(  ) const;
            real_function_type real_function_value( &::SireMaths::Complex::real );
            
            Complex_exposer.def( 
                "real"
                , real_function_value );
        
        }
        { //::SireMaths::Complex::rect
        
            typedef ::SireMaths::Complex ( *rect_function_type )( double,double );
            rect_function_type rect_function_value( &::SireMaths::Complex::rect );
            
            Complex_exposer.def( 
                "rect"
                , rect_function_value
                , ( bp::arg("x"), bp::arg("y") ) );
        
        }
        { //::SireMaths::Complex::setImag
        
            typedef void ( ::SireMaths::Complex::*setImag_function_type )( double ) ;
            setImag_function_type setImag_function_value( &::SireMaths::Complex::setImag );
            
            Complex_exposer.def( 
                "setImag"
                , setImag_function_value
                , ( bp::arg("y") ) );
        
        }
        { //::SireMaths::Complex::setPolar
        
            typedef void ( ::SireMaths::Complex::*setPolar_function_type )( double,double ) ;
            setPolar_function_type setPolar_function_value( &::SireMaths::Complex::setPolar );
            
            Complex_exposer.def( 
                "setPolar"
                , setPolar_function_value
                , ( bp::arg("r"), bp::arg("theta") ) );
        
        }
        { //::SireMaths::Complex::setReal
        
            typedef void ( ::SireMaths::Complex::*setReal_function_type )( double ) ;
            setReal_function_type setReal_function_value( &::SireMaths::Complex::setReal );
            
            Complex_exposer.def( 
                "setReal"
                , setReal_function_value
                , ( bp::arg("x") ) );
        
        }
        { //::SireMaths::Complex::setRectangular
        
            typedef void ( ::SireMaths::Complex::*setRectangular_function_type )( double,double ) ;
            setRectangular_function_type setRectangular_function_value( &::SireMaths::Complex::setRectangular );
            
            Complex_exposer.def( 
                "setRectangular"
                , setRectangular_function_value
                , ( bp::arg("x"), bp::arg("y") ) );
        
        }
        { //::SireMaths::Complex::toString
        
            typedef ::QString ( ::SireMaths::Complex::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMaths::Complex::toString );
            
            Complex_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMaths::Complex::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMaths::Complex::typeName );
            
            Complex_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMaths::Complex::what
        
            typedef char const * ( ::SireMaths::Complex::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireMaths::Complex::what );
            
            Complex_exposer.def( 
                "what"
                , what_function_value );
        
        }
        Complex_exposer.staticmethod( "polar" );
        Complex_exposer.staticmethod( "rect" );
        Complex_exposer.staticmethod( "typeName" );
        Complex_exposer.def( bp::other< double >() * bp::self );
        Complex_exposer.def( bp::self * bp::other< double >() );
        Complex_exposer.def( bp::self * bp::self );
        Complex_exposer.def( bp::other< double >() + bp::self );
        Complex_exposer.def( bp::self + bp::other< double >() );
        Complex_exposer.def( bp::self + bp::self );
        Complex_exposer.def( bp::other< double >() - bp::self );
        Complex_exposer.def( bp::self - bp::other< double >() );
        Complex_exposer.def( bp::self - bp::self );
        Complex_exposer.def( bp::other< double >() / bp::self );
        Complex_exposer.def( bp::self / bp::other< double >() );
        Complex_exposer.def( bp::self / bp::self );
        Complex_exposer.def( "__copy__", &__copy__);
        Complex_exposer.def( "__deepcopy__", &__copy__);
        Complex_exposer.def( "clone", &__copy__);
        Complex_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMaths::Complex >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Complex_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMaths::Complex >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Complex_exposer.def( "__str__", &__str__< ::SireMaths::Complex > );
        Complex_exposer.def( "__repr__", &__str__< ::SireMaths::Complex > );
    }

}