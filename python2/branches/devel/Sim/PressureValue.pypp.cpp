// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "PressureValue.pypp.hpp"

namespace bp = boost::python;

#include "dimensions.h"

#include "dimensions.hpp"

#include <cmath>

#include "dimensions.h"

SireSim::DimensionValue<SireSim::detail::PressureData> __copy__(const SireSim::DimensionValue<SireSim::detail::PressureData> &other){ return SireSim::DimensionValue<SireSim::detail::PressureData>(other); }

#include "Helpers/str.hpp"

void register_PressureValue_class(){

    { //::SireSim::DimensionValue< SireSim::detail::PressureData >
        typedef bp::class_< SireSim::DimensionValue< SireSim::detail::PressureData >, bp::bases< SireSim::Value > > PressureValue_exposer_t;
        PressureValue_exposer_t PressureValue_exposer = PressureValue_exposer_t( "PressureValue", bp::init< >() );
        bp::scope PressureValue_scope( PressureValue_exposer );
        PressureValue_exposer.def( bp::init< double, bp::optional< QString > >(( bp::arg("value"), bp::arg("units")=::QString( ) )) );
        PressureValue_exposer.def( bp::init< QString >(( bp::arg("value") )) );
        PressureValue_exposer.def( bp::init< QDomElement >(( bp::arg("elem") )) );
        PressureValue_exposer.def( bp::init< SireSim::DimensionValue< SireSim::detail::PressureData > const & >(( bp::arg("other") )) );
        PressureValue_exposer.def( bp::self != bp::self );
        { //::SireSim::DimensionValue< SireSim::detail::PressureData >::operator=
        
            typedef SireSim::DimensionValue< SireSim::detail::PressureData > exported_class_t;
            typedef ::SireSim::DimensionValue< SireSim::detail::PressureData > & ( ::SireSim::DimensionValue< SireSim::detail::PressureData >::*assign_function_type )( ::SireSim::DimensionValue< SireSim::detail::PressureData > const & ) ;
            assign_function_type assign_function_value( &::SireSim::DimensionValue< SireSim::detail::PressureData >::operator= );
            
            PressureValue_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        PressureValue_exposer.def( bp::self == bp::self );
        { //::SireSim::DimensionValue< SireSim::detail::PressureData >::supportedUnits
        
            typedef SireSim::DimensionValue< SireSim::detail::PressureData > exported_class_t;
            typedef ::QStringList ( *supportedUnits_function_type )(  );
            supportedUnits_function_type supportedUnits_function_value( &::SireSim::DimensionValue< SireSim::detail::PressureData >::supportedUnits );
            
            PressureValue_exposer.def( 
                "supportedUnits"
                , supportedUnits_function_value );
        
        }
        { //::SireSim::DimensionValue< SireSim::detail::PressureData >::typeName
        
            typedef SireSim::DimensionValue< SireSim::detail::PressureData > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireSim::DimensionValue< SireSim::detail::PressureData >::typeName );
            
            PressureValue_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireSim::DimensionValue< SireSim::detail::PressureData >::value
        
            typedef SireSim::DimensionValue< SireSim::detail::PressureData > exported_class_t;
            typedef double ( ::SireSim::DimensionValue< SireSim::detail::PressureData >::*value_function_type )(  ) const;
            value_function_type value_function_value( &::SireSim::DimensionValue< SireSim::detail::PressureData >::value );
            
            PressureValue_exposer.def( 
                "value"
                , value_function_value );
        
        }
        { //::SireSim::DimensionValue< SireSim::detail::PressureData >::value
        
            typedef SireSim::DimensionValue< SireSim::detail::PressureData > exported_class_t;
            typedef double ( ::SireSim::DimensionValue< SireSim::detail::PressureData >::*value_function_type )( ::QString ) const;
            value_function_type value_function_value( &::SireSim::DimensionValue< SireSim::detail::PressureData >::value );
            
            PressureValue_exposer.def( 
                "value"
                , value_function_value
                , ( bp::arg("units") ) );
        
        }
        { //::SireSim::DimensionValue< SireSim::detail::PressureData >::valueWithMaximum
        
            typedef SireSim::DimensionValue< SireSim::detail::PressureData > exported_class_t;
            typedef ::SireSim::DimensionValue< SireSim::detail::PressureData > ( *valueWithMaximum_function_type )( double,double,::QString );
            valueWithMaximum_function_type valueWithMaximum_function_value( &::SireSim::DimensionValue< SireSim::detail::PressureData >::valueWithMaximum );
            
            PressureValue_exposer.def( 
                "valueWithMaximum"
                , valueWithMaximum_function_value
                , ( bp::arg("value"), bp::arg("maximum"), bp::arg("units")=::QString( ) ) );
        
        }
        { //::SireSim::DimensionValue< SireSim::detail::PressureData >::valueWithMinimum
        
            typedef SireSim::DimensionValue< SireSim::detail::PressureData > exported_class_t;
            typedef ::SireSim::DimensionValue< SireSim::detail::PressureData > ( *valueWithMinimum_function_type )( double,double,::QString );
            valueWithMinimum_function_type valueWithMinimum_function_value( &::SireSim::DimensionValue< SireSim::detail::PressureData >::valueWithMinimum );
            
            PressureValue_exposer.def( 
                "valueWithMinimum"
                , valueWithMinimum_function_value
                , ( bp::arg("value"), bp::arg("minimum"), bp::arg("units")=::QString( ) ) );
        
        }
        { //::SireSim::DimensionValue< SireSim::detail::PressureData >::valueWithRange
        
            typedef SireSim::DimensionValue< SireSim::detail::PressureData > exported_class_t;
            typedef ::SireSim::DimensionValue< SireSim::detail::PressureData > ( *valueWithRange_function_type )( double,double,double,::QString );
            valueWithRange_function_type valueWithRange_function_value( &::SireSim::DimensionValue< SireSim::detail::PressureData >::valueWithRange );
            
            PressureValue_exposer.def( 
                "valueWithRange"
                , valueWithRange_function_value
                , ( bp::arg("value"), bp::arg("minimum"), bp::arg("maximum"), bp::arg("units")=::QString( ) ) );
        
        }
        { //::SireSim::DimensionValue< SireSim::detail::PressureData >::what
        
            typedef SireSim::DimensionValue< SireSim::detail::PressureData > exported_class_t;
            typedef char const * ( ::SireSim::DimensionValue< SireSim::detail::PressureData >::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireSim::DimensionValue< SireSim::detail::PressureData >::what );
            
            PressureValue_exposer.def( 
                "what"
                , what_function_value );
        
        }
        PressureValue_exposer.staticmethod( "supportedUnits" );
        PressureValue_exposer.staticmethod( "typeName" );
        PressureValue_exposer.staticmethod( "valueWithMaximum" );
        PressureValue_exposer.staticmethod( "valueWithMinimum" );
        PressureValue_exposer.staticmethod( "valueWithRange" );
        PressureValue_exposer.def( "__copy__", &__copy__);
        PressureValue_exposer.def( "__deepcopy__", &__copy__);
        PressureValue_exposer.def( "clone", &__copy__);
        PressureValue_exposer.def( "__str__", &__str__< ::SireSim::DimensionValue<SireSim::detail::PressureData> > );
        PressureValue_exposer.def( "__repr__", &__str__< ::SireSim::DimensionValue<SireSim::detail::PressureData> > );
    }

}