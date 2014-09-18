// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ChainStringProperty.pypp.hpp"

namespace bp = boost::python;

#include "chainproperty.hpp"

#include "chainproperty.hpp"

#include "SireMaths/vector.h"

SireMol::ChainProperty<QString> __copy__(const SireMol::ChainProperty<QString> &other){ return SireMol::ChainProperty<QString>(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_ChainStringProperty_class(){

    { //::SireMol::ChainProperty< QString >
        typedef bp::class_< SireMol::ChainProperty< QString >, bp::bases< SireMol::ChainProp, SireMol::MolViewProperty, SireBase::Property > > ChainStringProperty_exposer_t;
        ChainStringProperty_exposer_t ChainStringProperty_exposer = ChainStringProperty_exposer_t( "ChainStringProperty", bp::init< >() );
        bp::scope ChainStringProperty_scope( ChainStringProperty_exposer );
        ChainStringProperty_exposer.def( bp::init< SireMol::MoleculeInfoData const & >(( bp::arg("molinfo") )) );
        ChainStringProperty_exposer.def( bp::init< QVector< QString > const & >(( bp::arg("values") )) );
        ChainStringProperty_exposer.def( bp::init< SireMol::ChainProperty< QString > const & >(( bp::arg("other") )) );
        { //::SireMol::ChainProperty< QString >::array
        
            typedef SireMol::ChainProperty< QString > exported_class_t;
            typedef ::QVector< QString > const & ( ::SireMol::ChainProperty< QString >::*array_function_type )(  ) const;
            array_function_type array_function_value( &::SireMol::ChainProperty< QString >::array );
            
            ChainStringProperty_exposer.def( 
                "array"
                , array_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::ChainProperty< QString >::assertCanConvert
        
            typedef SireMol::ChainProperty< QString > exported_class_t;
            typedef void ( ::SireMol::ChainProperty< QString >::*assertCanConvert_function_type )( ::QVariant const & ) const;
            assertCanConvert_function_type assertCanConvert_function_value( &::SireMol::ChainProperty< QString >::assertCanConvert );
            
            ChainStringProperty_exposer.def( 
                "assertCanConvert"
                , assertCanConvert_function_value
                , ( bp::arg("value") ) );
        
        }
        { //::SireMol::ChainProperty< QString >::assignFrom
        
            typedef SireMol::ChainProperty< QString > exported_class_t;
            typedef void ( ::SireMol::ChainProperty< QString >::*assignFrom_function_type )( ::SireMol::ChainProperty< QVariant > const & ) ;
            assignFrom_function_type assignFrom_function_value( &::SireMol::ChainProperty< QString >::assignFrom );
            
            ChainStringProperty_exposer.def( 
                "assignFrom"
                , assignFrom_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireMol::ChainProperty< QString >::at
        
            typedef SireMol::ChainProperty< QString > exported_class_t;
            typedef ::QString const & ( ::SireMol::ChainProperty< QString >::*at_function_type )( ::SireMol::ChainIdx const & ) const;
            at_function_type at_function_value( &::SireMol::ChainProperty< QString >::at );
            
            ChainStringProperty_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("chainidx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::ChainProperty< QString >::canConvert
        
            typedef SireMol::ChainProperty< QString > exported_class_t;
            typedef bool ( ::SireMol::ChainProperty< QString >::*canConvert_function_type )( ::QVariant const & ) const;
            canConvert_function_type canConvert_function_value( &::SireMol::ChainProperty< QString >::canConvert );
            
            ChainStringProperty_exposer.def( 
                "canConvert"
                , canConvert_function_value
                , ( bp::arg("value") ) );
        
        }
        { //::SireMol::ChainProperty< QString >::count
        
            typedef SireMol::ChainProperty< QString > exported_class_t;
            typedef int ( ::SireMol::ChainProperty< QString >::*count_function_type )(  ) const;
            count_function_type count_function_value( &::SireMol::ChainProperty< QString >::count );
            
            ChainStringProperty_exposer.def( 
                "count"
                , count_function_value );
        
        }
        { //::SireMol::ChainProperty< QString >::fromVariant
        
            typedef SireMol::ChainProperty< QString > exported_class_t;
            typedef ::SireMol::ChainProperty< QString > ( *fromVariant_function_type )( ::SireMol::ChainProperty< QVariant > const & );
            fromVariant_function_type fromVariant_function_value( &::SireMol::ChainProperty< QString >::fromVariant );
            
            ChainStringProperty_exposer.def( 
                "fromVariant"
                , fromVariant_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireMol::ChainProperty< QString >::get
        
            typedef SireMol::ChainProperty< QString > exported_class_t;
            typedef ::QString const & ( ::SireMol::ChainProperty< QString >::*get_function_type )( ::SireMol::ChainIdx const & ) const;
            get_function_type get_function_value( &::SireMol::ChainProperty< QString >::get );
            
            ChainStringProperty_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("chainidx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::ChainProperty< QString >::isCompatibleWith
        
            typedef SireMol::ChainProperty< QString > exported_class_t;
            typedef bool ( ::SireMol::ChainProperty< QString >::*isCompatibleWith_function_type )( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::ChainProperty< QString >::isCompatibleWith );
            
            ChainStringProperty_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") ) );
        
        }
        { //::SireMol::ChainProperty< QString >::isEmpty
        
            typedef SireMol::ChainProperty< QString > exported_class_t;
            typedef bool ( ::SireMol::ChainProperty< QString >::*isEmpty_function_type )(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::ChainProperty< QString >::isEmpty );
            
            ChainStringProperty_exposer.def( 
                "isEmpty"
                , isEmpty_function_value );
        
        }
        { //::SireMol::ChainProperty< QString >::nChains
        
            typedef SireMol::ChainProperty< QString > exported_class_t;
            typedef int ( ::SireMol::ChainProperty< QString >::*nChains_function_type )(  ) const;
            nChains_function_type nChains_function_value( &::SireMol::ChainProperty< QString >::nChains );
            
            ChainStringProperty_exposer.def( 
                "nChains"
                , nChains_function_value );
        
        }
        ChainStringProperty_exposer.def( bp::self != bp::self );
        { //::SireMol::ChainProperty< QString >::operator=
        
            typedef SireMol::ChainProperty< QString > exported_class_t;
            typedef ::SireMol::ChainProperty< QString > & ( ::SireMol::ChainProperty< QString >::*assign_function_type )( ::SireMol::ChainProperty< QString > const & ) ;
            assign_function_type assign_function_value( &::SireMol::ChainProperty< QString >::operator= );
            
            ChainStringProperty_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        ChainStringProperty_exposer.def( bp::self == bp::self );
        { //::SireMol::ChainProperty< QString >::operator[]
        
            typedef SireMol::ChainProperty< QString > exported_class_t;
            typedef ::QString const & ( ::SireMol::ChainProperty< QString >::*__getitem___function_type )( ::SireMol::ChainIdx const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::ChainProperty< QString >::operator[] );
            
            ChainStringProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("chainidx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::ChainProperty< QString >::set
        
            typedef SireMol::ChainProperty< QString > exported_class_t;
            typedef ::SireMol::ChainProperty< QString > & ( ::SireMol::ChainProperty< QString >::*set_function_type )( ::SireMol::ChainIdx,::QString const & ) ;
            set_function_type set_function_value( &::SireMol::ChainProperty< QString >::set );
            
            ChainStringProperty_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("chainidx"), bp::arg("value") )
                , bp::return_self< >() );
        
        }
        { //::SireMol::ChainProperty< QString >::size
        
            typedef SireMol::ChainProperty< QString > exported_class_t;
            typedef int ( ::SireMol::ChainProperty< QString >::*size_function_type )(  ) const;
            size_function_type size_function_value( &::SireMol::ChainProperty< QString >::size );
            
            ChainStringProperty_exposer.def( 
                "size"
                , size_function_value );
        
        }
        { //::SireMol::ChainProperty< QString >::toString
        
            typedef SireMol::ChainProperty< QString > exported_class_t;
            typedef ::QString ( ::SireMol::ChainProperty< QString >::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMol::ChainProperty< QString >::toString );
            
            ChainStringProperty_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMol::ChainProperty< QString >::toVariant
        
            typedef SireMol::ChainProperty< QString > exported_class_t;
            typedef ::SireMol::ChainProperty< QVariant > ( ::SireMol::ChainProperty< QString >::*toVariant_function_type )(  ) const;
            toVariant_function_type toVariant_function_value( &::SireMol::ChainProperty< QString >::toVariant );
            
            ChainStringProperty_exposer.def( 
                "toVariant"
                , toVariant_function_value );
        
        }
        { //::SireMol::ChainProperty< QString >::typeName
        
            typedef SireMol::ChainProperty< QString > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::ChainProperty< QString >::typeName );
            
            ChainStringProperty_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        ChainStringProperty_exposer.staticmethod( "fromVariant" );
        ChainStringProperty_exposer.staticmethod( "typeName" );
        ChainStringProperty_exposer.def( "__copy__", &__copy__);
        ChainStringProperty_exposer.def( "__deepcopy__", &__copy__);
        ChainStringProperty_exposer.def( "clone", &__copy__);
        ChainStringProperty_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::ChainProperty<QString> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ChainStringProperty_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::ChainProperty<QString> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ChainStringProperty_exposer.def( "__str__", &__str__< ::SireMol::ChainProperty<QString> > );
        ChainStringProperty_exposer.def( "__repr__", &__str__< ::SireMol::ChainProperty<QString> > );
        ChainStringProperty_exposer.def( "__len__", &__len_size< ::SireMol::ChainProperty<QString> > );
    }

}
