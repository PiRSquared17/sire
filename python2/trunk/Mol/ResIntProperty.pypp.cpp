// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ResIntProperty.pypp.hpp"

namespace bp = boost::python;

#include "resproperty.hpp"

#include "resproperty.hpp"

#include "SireMaths/vector.h"

SireMol::ResProperty<long long> __copy__(const SireMol::ResProperty<long long> &other){ return SireMol::ResProperty<long long>(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_ResIntProperty_class(){

    { //::SireMol::ResProperty< long long >
        typedef bp::class_< SireMol::ResProperty< long long >, bp::bases< SireMol::ResProp, SireMol::MolViewProperty, SireBase::Property > > ResIntProperty_exposer_t;
        ResIntProperty_exposer_t ResIntProperty_exposer = ResIntProperty_exposer_t( "ResIntProperty", bp::init< >() );
        bp::scope ResIntProperty_scope( ResIntProperty_exposer );
        ResIntProperty_exposer.def( bp::init< SireMol::MoleculeInfoData const & >(( bp::arg("molinfo") )) );
        ResIntProperty_exposer.def( bp::init< QVector< long long > const & >(( bp::arg("values") )) );
        ResIntProperty_exposer.def( bp::init< SireMol::ResProperty< long long > const & >(( bp::arg("other") )) );
        { //::SireMol::ResProperty< long long >::array
        
            typedef SireMol::ResProperty< long long > exported_class_t;
            typedef ::QVector< long long > const & ( ::SireMol::ResProperty< long long >::*array_function_type )(  ) const;
            array_function_type array_function_value( &::SireMol::ResProperty< long long >::array );
            
            ResIntProperty_exposer.def( 
                "array"
                , array_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::ResProperty< long long >::assertCanConvert
        
            typedef SireMol::ResProperty< long long > exported_class_t;
            typedef void ( ::SireMol::ResProperty< long long >::*assertCanConvert_function_type )( ::QVariant const & ) const;
            assertCanConvert_function_type assertCanConvert_function_value( &::SireMol::ResProperty< long long >::assertCanConvert );
            
            ResIntProperty_exposer.def( 
                "assertCanConvert"
                , assertCanConvert_function_value
                , ( bp::arg("value") ) );
        
        }
        { //::SireMol::ResProperty< long long >::assignFrom
        
            typedef SireMol::ResProperty< long long > exported_class_t;
            typedef void ( ::SireMol::ResProperty< long long >::*assignFrom_function_type )( ::SireMol::ResProperty< QVariant > const & ) ;
            assignFrom_function_type assignFrom_function_value( &::SireMol::ResProperty< long long >::assignFrom );
            
            ResIntProperty_exposer.def( 
                "assignFrom"
                , assignFrom_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireMol::ResProperty< long long >::at
        
            typedef SireMol::ResProperty< long long > exported_class_t;
            typedef long long int const & ( ::SireMol::ResProperty< long long >::*at_function_type )( ::SireMol::ResIdx const & ) const;
            at_function_type at_function_value( &::SireMol::ResProperty< long long >::at );
            
            ResIntProperty_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("residx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::ResProperty< long long >::canConvert
        
            typedef SireMol::ResProperty< long long > exported_class_t;
            typedef bool ( ::SireMol::ResProperty< long long >::*canConvert_function_type )( ::QVariant const & ) const;
            canConvert_function_type canConvert_function_value( &::SireMol::ResProperty< long long >::canConvert );
            
            ResIntProperty_exposer.def( 
                "canConvert"
                , canConvert_function_value
                , ( bp::arg("value") ) );
        
        }
        { //::SireMol::ResProperty< long long >::count
        
            typedef SireMol::ResProperty< long long > exported_class_t;
            typedef int ( ::SireMol::ResProperty< long long >::*count_function_type )(  ) const;
            count_function_type count_function_value( &::SireMol::ResProperty< long long >::count );
            
            ResIntProperty_exposer.def( 
                "count"
                , count_function_value );
        
        }
        { //::SireMol::ResProperty< long long >::fromVariant
        
            typedef SireMol::ResProperty< long long > exported_class_t;
            typedef ::SireMol::ResProperty< long long > ( *fromVariant_function_type )( ::SireMol::ResProperty< QVariant > const & );
            fromVariant_function_type fromVariant_function_value( &::SireMol::ResProperty< long long >::fromVariant );
            
            ResIntProperty_exposer.def( 
                "fromVariant"
                , fromVariant_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireMol::ResProperty< long long >::get
        
            typedef SireMol::ResProperty< long long > exported_class_t;
            typedef long long int const & ( ::SireMol::ResProperty< long long >::*get_function_type )( ::SireMol::ResIdx const & ) const;
            get_function_type get_function_value( &::SireMol::ResProperty< long long >::get );
            
            ResIntProperty_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("residx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::ResProperty< long long >::isCompatibleWith
        
            typedef SireMol::ResProperty< long long > exported_class_t;
            typedef bool ( ::SireMol::ResProperty< long long >::*isCompatibleWith_function_type )( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::ResProperty< long long >::isCompatibleWith );
            
            ResIntProperty_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") ) );
        
        }
        { //::SireMol::ResProperty< long long >::isEmpty
        
            typedef SireMol::ResProperty< long long > exported_class_t;
            typedef bool ( ::SireMol::ResProperty< long long >::*isEmpty_function_type )(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::ResProperty< long long >::isEmpty );
            
            ResIntProperty_exposer.def( 
                "isEmpty"
                , isEmpty_function_value );
        
        }
        { //::SireMol::ResProperty< long long >::nResidues
        
            typedef SireMol::ResProperty< long long > exported_class_t;
            typedef int ( ::SireMol::ResProperty< long long >::*nResidues_function_type )(  ) const;
            nResidues_function_type nResidues_function_value( &::SireMol::ResProperty< long long >::nResidues );
            
            ResIntProperty_exposer.def( 
                "nResidues"
                , nResidues_function_value );
        
        }
        ResIntProperty_exposer.def( bp::self != bp::self );
        { //::SireMol::ResProperty< long long >::operator=
        
            typedef SireMol::ResProperty< long long > exported_class_t;
            typedef ::SireMol::ResProperty< long long > & ( ::SireMol::ResProperty< long long >::*assign_function_type )( ::SireMol::ResProperty< long long > const & ) ;
            assign_function_type assign_function_value( &::SireMol::ResProperty< long long >::operator= );
            
            ResIntProperty_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        ResIntProperty_exposer.def( bp::self == bp::self );
        { //::SireMol::ResProperty< long long >::operator[]
        
            typedef SireMol::ResProperty< long long > exported_class_t;
            typedef long long int const & ( ::SireMol::ResProperty< long long >::*__getitem___function_type )( ::SireMol::ResIdx const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::ResProperty< long long >::operator[] );
            
            ResIntProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("residx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::ResProperty< long long >::set
        
            typedef SireMol::ResProperty< long long > exported_class_t;
            typedef ::SireMol::ResProperty< long long > & ( ::SireMol::ResProperty< long long >::*set_function_type )( ::SireMol::ResIdx,long long int const & ) ;
            set_function_type set_function_value( &::SireMol::ResProperty< long long >::set );
            
            ResIntProperty_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("residx"), bp::arg("value") )
                , bp::return_self< >() );
        
        }
        { //::SireMol::ResProperty< long long >::size
        
            typedef SireMol::ResProperty< long long > exported_class_t;
            typedef int ( ::SireMol::ResProperty< long long >::*size_function_type )(  ) const;
            size_function_type size_function_value( &::SireMol::ResProperty< long long >::size );
            
            ResIntProperty_exposer.def( 
                "size"
                , size_function_value );
        
        }
        { //::SireMol::ResProperty< long long >::toString
        
            typedef SireMol::ResProperty< long long > exported_class_t;
            typedef ::QString ( ::SireMol::ResProperty< long long >::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMol::ResProperty< long long >::toString );
            
            ResIntProperty_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMol::ResProperty< long long >::toVariant
        
            typedef SireMol::ResProperty< long long > exported_class_t;
            typedef ::SireMol::ResProperty< QVariant > ( ::SireMol::ResProperty< long long >::*toVariant_function_type )(  ) const;
            toVariant_function_type toVariant_function_value( &::SireMol::ResProperty< long long >::toVariant );
            
            ResIntProperty_exposer.def( 
                "toVariant"
                , toVariant_function_value );
        
        }
        { //::SireMol::ResProperty< long long >::typeName
        
            typedef SireMol::ResProperty< long long > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::ResProperty< long long >::typeName );
            
            ResIntProperty_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        ResIntProperty_exposer.staticmethod( "fromVariant" );
        ResIntProperty_exposer.staticmethod( "typeName" );
        ResIntProperty_exposer.def( "__copy__", &__copy__);
        ResIntProperty_exposer.def( "__deepcopy__", &__copy__);
        ResIntProperty_exposer.def( "clone", &__copy__);
        ResIntProperty_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::ResProperty<long long> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ResIntProperty_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::ResProperty<long long> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ResIntProperty_exposer.def( "__str__", &__str__< ::SireMol::ResProperty<long long> > );
        ResIntProperty_exposer.def( "__repr__", &__str__< ::SireMol::ResProperty<long long> > );
    }

}
