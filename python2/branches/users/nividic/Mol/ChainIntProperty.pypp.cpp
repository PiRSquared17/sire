// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ChainIntProperty.pypp.hpp"

namespace bp = boost::python;

#include "chainproperty.hpp"

#include "chainproperty.hpp"

#include "SireMaths/vector.h"

SireMol::ChainProperty<long long> __copy__(const SireMol::ChainProperty<long long> &other){ return SireMol::ChainProperty<long long>(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_ChainIntProperty_class(){

    { //::SireMol::ChainProperty< long long int >
        typedef bp::class_< SireMol::ChainProperty< long long int >, bp::bases< SireMol::ChainProp, SireMol::MolViewProperty, SireBase::Property > > ChainIntProperty_exposer_t;
        ChainIntProperty_exposer_t ChainIntProperty_exposer = ChainIntProperty_exposer_t( "ChainIntProperty" );
        bp::scope ChainIntProperty_scope( ChainIntProperty_exposer );
        ChainIntProperty_exposer.def( bp::init< >() );
        ChainIntProperty_exposer.def( bp::init< SireMol::MoleculeInfoData const & >(( bp::arg("molinfo") )) );
        ChainIntProperty_exposer.def( bp::init< QVector< long long int > const & >(( bp::arg("values") )) );
        ChainIntProperty_exposer.def( bp::init< SireMol::ChainProperty< long long int > const & >(( bp::arg("other") )) );
        { //::SireMol::ChainProperty< long long int >::array
        
            typedef SireMol::ChainProperty< long long int > exported_class_t;
            typedef ::QVector< long long int > const & ( ::SireMol::ChainProperty< long long int >::*array_function_type )(  ) const;
            array_function_type array_function_value( &::SireMol::ChainProperty< long long int >::array );
            
            ChainIntProperty_exposer.def( 
                "array"
                , array_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::ChainProperty< long long int >::assertCanConvert
        
            typedef SireMol::ChainProperty< long long int > exported_class_t;
            typedef void ( ::SireMol::ChainProperty< long long int >::*assertCanConvert_function_type )( ::QVariant const & ) const;
            assertCanConvert_function_type assertCanConvert_function_value( &::SireMol::ChainProperty< long long int >::assertCanConvert );
            
            ChainIntProperty_exposer.def( 
                "assertCanConvert"
                , assertCanConvert_function_value
                , ( bp::arg("value") ) );
        
        }
        { //::SireMol::ChainProperty< long long int >::assignFrom
        
            typedef SireMol::ChainProperty< long long int > exported_class_t;
            typedef void ( ::SireMol::ChainProperty< long long int >::*assignFrom_function_type )( ::SireMol::ChainProperty< QVariant > const & ) ;
            assignFrom_function_type assignFrom_function_value( &::SireMol::ChainProperty< long long int >::assignFrom );
            
            ChainIntProperty_exposer.def( 
                "assignFrom"
                , assignFrom_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireMol::ChainProperty< long long int >::at
        
            typedef SireMol::ChainProperty< long long int > exported_class_t;
            typedef long long int const & ( ::SireMol::ChainProperty< long long int >::*at_function_type )( ::SireMol::ChainIdx const & ) const;
            at_function_type at_function_value( &::SireMol::ChainProperty< long long int >::at );
            
            ChainIntProperty_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("chainidx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::ChainProperty< long long int >::canConvert
        
            typedef SireMol::ChainProperty< long long int > exported_class_t;
            typedef bool ( ::SireMol::ChainProperty< long long int >::*canConvert_function_type )( ::QVariant const & ) const;
            canConvert_function_type canConvert_function_value( &::SireMol::ChainProperty< long long int >::canConvert );
            
            ChainIntProperty_exposer.def( 
                "canConvert"
                , canConvert_function_value
                , ( bp::arg("value") ) );
        
        }
        { //::SireMol::ChainProperty< long long int >::count
        
            typedef SireMol::ChainProperty< long long int > exported_class_t;
            typedef int ( ::SireMol::ChainProperty< long long int >::*count_function_type )(  ) const;
            count_function_type count_function_value( &::SireMol::ChainProperty< long long int >::count );
            
            ChainIntProperty_exposer.def( 
                "count"
                , count_function_value );
        
        }
        { //::SireMol::ChainProperty< long long int >::fromVariant
        
            typedef SireMol::ChainProperty< long long int > exported_class_t;
            typedef ::SireMol::ChainProperty< long long int > ( *fromVariant_function_type )( ::SireMol::ChainProperty< QVariant > const & );
            fromVariant_function_type fromVariant_function_value( &::SireMol::ChainProperty< long long int >::fromVariant );
            
            ChainIntProperty_exposer.def( 
                "fromVariant"
                , fromVariant_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireMol::ChainProperty< long long int >::get
        
            typedef SireMol::ChainProperty< long long int > exported_class_t;
            typedef long long int const & ( ::SireMol::ChainProperty< long long int >::*get_function_type )( ::SireMol::ChainIdx const & ) const;
            get_function_type get_function_value( &::SireMol::ChainProperty< long long int >::get );
            
            ChainIntProperty_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("chainidx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::ChainProperty< long long int >::isCompatibleWith
        
            typedef SireMol::ChainProperty< long long int > exported_class_t;
            typedef bool ( ::SireMol::ChainProperty< long long int >::*isCompatibleWith_function_type )( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::ChainProperty< long long int >::isCompatibleWith );
            
            ChainIntProperty_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") ) );
        
        }
        { //::SireMol::ChainProperty< long long int >::isEmpty
        
            typedef SireMol::ChainProperty< long long int > exported_class_t;
            typedef bool ( ::SireMol::ChainProperty< long long int >::*isEmpty_function_type )(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::ChainProperty< long long int >::isEmpty );
            
            ChainIntProperty_exposer.def( 
                "isEmpty"
                , isEmpty_function_value );
        
        }
        { //::SireMol::ChainProperty< long long int >::nChains
        
            typedef SireMol::ChainProperty< long long int > exported_class_t;
            typedef int ( ::SireMol::ChainProperty< long long int >::*nChains_function_type )(  ) const;
            nChains_function_type nChains_function_value( &::SireMol::ChainProperty< long long int >::nChains );
            
            ChainIntProperty_exposer.def( 
                "nChains"
                , nChains_function_value );
        
        }
        ChainIntProperty_exposer.def( bp::self != bp::self );
        { //::SireMol::ChainProperty< long long int >::operator=
        
            typedef SireMol::ChainProperty< long long int > exported_class_t;
            typedef ::SireMol::ChainProperty< long long int > & ( ::SireMol::ChainProperty< long long int >::*assign_function_type )( ::SireMol::ChainProperty< long long int > const & ) ;
            assign_function_type assign_function_value( &::SireMol::ChainProperty< long long int >::operator= );
            
            ChainIntProperty_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        ChainIntProperty_exposer.def( bp::self == bp::self );
        { //::SireMol::ChainProperty< long long int >::operator[]
        
            typedef SireMol::ChainProperty< long long int > exported_class_t;
            typedef long long int const & ( ::SireMol::ChainProperty< long long int >::*__getitem___function_type )( ::SireMol::ChainIdx const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::ChainProperty< long long int >::operator[] );
            
            ChainIntProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("chainidx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::ChainProperty< long long int >::set
        
            typedef SireMol::ChainProperty< long long int > exported_class_t;
            typedef ::SireMol::ChainProperty< long long int > & ( ::SireMol::ChainProperty< long long int >::*set_function_type )( ::SireMol::ChainIdx,long long int const & ) ;
            set_function_type set_function_value( &::SireMol::ChainProperty< long long int >::set );
            
            ChainIntProperty_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("chainidx"), bp::arg("value") )
                , bp::return_self< >() );
        
        }
        { //::SireMol::ChainProperty< long long int >::size
        
            typedef SireMol::ChainProperty< long long int > exported_class_t;
            typedef int ( ::SireMol::ChainProperty< long long int >::*size_function_type )(  ) const;
            size_function_type size_function_value( &::SireMol::ChainProperty< long long int >::size );
            
            ChainIntProperty_exposer.def( 
                "size"
                , size_function_value );
        
        }
        { //::SireMol::ChainProperty< long long int >::toString
        
            typedef SireMol::ChainProperty< long long int > exported_class_t;
            typedef ::QString ( ::SireMol::ChainProperty< long long int >::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMol::ChainProperty< long long int >::toString );
            
            ChainIntProperty_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMol::ChainProperty< long long int >::toVariant
        
            typedef SireMol::ChainProperty< long long int > exported_class_t;
            typedef ::SireMol::ChainProperty< QVariant > ( ::SireMol::ChainProperty< long long int >::*toVariant_function_type )(  ) const;
            toVariant_function_type toVariant_function_value( &::SireMol::ChainProperty< long long int >::toVariant );
            
            ChainIntProperty_exposer.def( 
                "toVariant"
                , toVariant_function_value );
        
        }
        { //::SireMol::ChainProperty< long long int >::typeName
        
            typedef SireMol::ChainProperty< long long int > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::ChainProperty< long long int >::typeName );
            
            ChainIntProperty_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        ChainIntProperty_exposer.staticmethod( "fromVariant" );
        ChainIntProperty_exposer.staticmethod( "typeName" );
        ChainIntProperty_exposer.def( "__copy__", &__copy__);
        ChainIntProperty_exposer.def( "__deepcopy__", &__copy__);
        ChainIntProperty_exposer.def( "clone", &__copy__);
        ChainIntProperty_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::ChainProperty<long long> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ChainIntProperty_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::ChainProperty<long long> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ChainIntProperty_exposer.def( "__str__", &__str__< ::SireMol::ChainProperty<long long> > );
        ChainIntProperty_exposer.def( "__repr__", &__str__< ::SireMol::ChainProperty<long long> > );
    }

}
