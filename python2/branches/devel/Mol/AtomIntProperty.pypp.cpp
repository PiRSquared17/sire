// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "AtomIntProperty.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMaths/vector.h"

#include "atomcharges.h"

#include "atomelements.h"

#include "atomenergies.h"

#include "atomforces.h"

#include "atommasses.h"

#include "atompolarisabilities.h"

#include "atomproperty.hpp"

#include "atomvelocities.h"

#include "atomproperty.hpp"

#include "SireMaths/vector.h"

SireMol::AtomProperty<long long> __copy__(const SireMol::AtomProperty<long long> &other){ return SireMol::AtomProperty<long long>(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_AtomIntProperty_class(){

    { //::SireMol::AtomProperty< long long int >
        typedef bp::class_< SireMol::AtomProperty< long long int >, bp::bases< SireMol::AtomProp, SireMol::MolViewProperty, SireBase::Property > > AtomIntProperty_exposer_t;
        AtomIntProperty_exposer_t AtomIntProperty_exposer = AtomIntProperty_exposer_t( "AtomIntProperty" );
        bp::scope AtomIntProperty_scope( AtomIntProperty_exposer );
        AtomIntProperty_exposer.def( bp::init< >() );
        AtomIntProperty_exposer.def( bp::init< SireMol::MoleculeInfoData const & >(( bp::arg("molinfo") )) );
        AtomIntProperty_exposer.def( bp::init< long long int const & >(( bp::arg("value") )) );
        AtomIntProperty_exposer.def( bp::init< SireBase::PackedArray2D< long long int > const & >(( bp::arg("values") )) );
        AtomIntProperty_exposer.def( bp::init< SireMol::AtomProperty< long long int > const & >(( bp::arg("other") )) );
        { //::SireMol::AtomProperty< long long int >::array
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef ::SireBase::PackedArray2D< long long int > const & ( ::SireMol::AtomProperty< long long int >::*array_function_type )(  ) const;
            array_function_type array_function_value( &::SireMol::AtomProperty< long long int >::array );
            
            AtomIntProperty_exposer.def( 
                "array"
                , array_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::AtomProperty< long long int >::assertCanConvert
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< long long int >::*assertCanConvert_function_type )( ::QVariant const & ) const;
            assertCanConvert_function_type assertCanConvert_function_value( &::SireMol::AtomProperty< long long int >::assertCanConvert );
            
            AtomIntProperty_exposer.def( 
                "assertCanConvert"
                , assertCanConvert_function_value
                , ( bp::arg("value") ) );
        
        }
        { //::SireMol::AtomProperty< long long int >::assignFrom
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< long long int >::*assignFrom_function_type )( ::SireMol::AtomProperty< QVariant > const & ) ;
            assignFrom_function_type assignFrom_function_value( &::SireMol::AtomProperty< long long int >::assignFrom );
            
            AtomIntProperty_exposer.def( 
                "assignFrom"
                , assignFrom_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireMol::AtomProperty< long long int >::at
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef ::SireBase::detail::PackedArray2D_Array< long long int > const & ( ::SireMol::AtomProperty< long long int >::*at_function_type )( ::SireMol::CGIdx ) const;
            at_function_type at_function_value( &::SireMol::AtomProperty< long long int >::at );
            
            AtomIntProperty_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("cgidx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::AtomProperty< long long int >::at
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef long long int const & ( ::SireMol::AtomProperty< long long int >::*at_function_type )( ::SireMol::CGAtomIdx const & ) const;
            at_function_type at_function_value( &::SireMol::AtomProperty< long long int >::at );
            
            AtomIntProperty_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("cgatomidx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::AtomProperty< long long int >::canConvert
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< long long int >::*canConvert_function_type )( ::QVariant const & ) const;
            canConvert_function_type canConvert_function_value( &::SireMol::AtomProperty< long long int >::canConvert );
            
            AtomIntProperty_exposer.def( 
                "canConvert"
                , canConvert_function_value
                , ( bp::arg("value") ) );
        
        }
        { //::SireMol::AtomProperty< long long int >::copyFrom
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< long long int >::*copyFrom_function_type )( ::QVector< long long int > const & ) ;
            copyFrom_function_type copyFrom_function_value( &::SireMol::AtomProperty< long long int >::copyFrom );
            
            AtomIntProperty_exposer.def( 
                "copyFrom"
                , copyFrom_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireMol::AtomProperty< long long int >::copyFrom
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< long long int >::*copyFrom_function_type )( ::QVector< long long int > const &,::SireMol::AtomSelection const & ) ;
            copyFrom_function_type copyFrom_function_value( &::SireMol::AtomProperty< long long int >::copyFrom );
            
            AtomIntProperty_exposer.def( 
                "copyFrom"
                , copyFrom_function_value
                , ( bp::arg("values"), bp::arg("selection") ) );
        
        }
        { //::SireMol::AtomProperty< long long int >::count
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< long long int >::*count_function_type )(  ) const;
            count_function_type count_function_value( &::SireMol::AtomProperty< long long int >::count );
            
            AtomIntProperty_exposer.def( 
                "count"
                , count_function_value );
        
        }
        { //::SireMol::AtomProperty< long long int >::fromVariant
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef ::SireMol::AtomProperty< long long int > ( *fromVariant_function_type )( ::SireMol::AtomProperty< QVariant > const & );
            fromVariant_function_type fromVariant_function_value( &::SireMol::AtomProperty< long long int >::fromVariant );
            
            AtomIntProperty_exposer.def( 
                "fromVariant"
                , fromVariant_function_value
                , ( bp::arg("variant") ) );
        
        }
        { //::SireMol::AtomProperty< long long int >::get
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef ::SireBase::detail::PackedArray2D_Array< long long int > const & ( ::SireMol::AtomProperty< long long int >::*get_function_type )( ::SireMol::CGIdx ) const;
            get_function_type get_function_value( &::SireMol::AtomProperty< long long int >::get );
            
            AtomIntProperty_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("cgidx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::AtomProperty< long long int >::get
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef long long int const & ( ::SireMol::AtomProperty< long long int >::*get_function_type )( ::SireMol::CGAtomIdx const & ) const;
            get_function_type get_function_value( &::SireMol::AtomProperty< long long int >::get );
            
            AtomIntProperty_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("cgatomidx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::AtomProperty< long long int >::isCompatibleWith
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< long long int >::*isCompatibleWith_function_type )( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::AtomProperty< long long int >::isCompatibleWith );
            
            AtomIntProperty_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") ) );
        
        }
        { //::SireMol::AtomProperty< long long int >::isEmpty
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< long long int >::*isEmpty_function_type )(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::AtomProperty< long long int >::isEmpty );
            
            AtomIntProperty_exposer.def( 
                "isEmpty"
                , isEmpty_function_value );
        
        }
        { //::SireMol::AtomProperty< long long int >::matchToSelection
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef ::SireMol::AtomProperty< long long int > ( ::SireMol::AtomProperty< long long int >::*matchToSelection_function_type )( ::SireMol::AtomSelection const & ) const;
            matchToSelection_function_type matchToSelection_function_value( &::SireMol::AtomProperty< long long int >::matchToSelection );
            
            AtomIntProperty_exposer.def( 
                "matchToSelection"
                , matchToSelection_function_value
                , ( bp::arg("selection") ) );
        
        }
        { //::SireMol::AtomProperty< long long int >::nAtoms
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< long long int >::*nAtoms_function_type )(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::AtomProperty< long long int >::nAtoms );
            
            AtomIntProperty_exposer.def( 
                "nAtoms"
                , nAtoms_function_value );
        
        }
        { //::SireMol::AtomProperty< long long int >::nAtoms
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< long long int >::*nAtoms_function_type )( ::SireMol::CGIdx ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::AtomProperty< long long int >::nAtoms );
            
            AtomIntProperty_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , ( bp::arg("cgidx") ) );
        
        }
        { //::SireMol::AtomProperty< long long int >::nCutGroups
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< long long int >::*nCutGroups_function_type )(  ) const;
            nCutGroups_function_type nCutGroups_function_value( &::SireMol::AtomProperty< long long int >::nCutGroups );
            
            AtomIntProperty_exposer.def( 
                "nCutGroups"
                , nCutGroups_function_value );
        
        }
        AtomIntProperty_exposer.def( bp::self != bp::self );
        { //::SireMol::AtomProperty< long long int >::operator=
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef ::SireMol::AtomProperty< long long int > & ( ::SireMol::AtomProperty< long long int >::*assign_function_type )( ::SireMol::AtomProperty< long long int > const & ) ;
            assign_function_type assign_function_value( &::SireMol::AtomProperty< long long int >::operator= );
            
            AtomIntProperty_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        AtomIntProperty_exposer.def( bp::self == bp::self );
        { //::SireMol::AtomProperty< long long int >::operator[]
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef ::SireBase::detail::PackedArray2D_Array< long long int > const & ( ::SireMol::AtomProperty< long long int >::*__getitem___function_type )( ::SireMol::CGIdx ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::AtomProperty< long long int >::operator[] );
            
            AtomIntProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("cgidx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::AtomProperty< long long int >::operator[]
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef long long int const & ( ::SireMol::AtomProperty< long long int >::*__getitem___function_type )( ::SireMol::CGAtomIdx const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::AtomProperty< long long int >::operator[] );
            
            AtomIntProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("cgatomidx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::AtomProperty< long long int >::set
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef ::SireMol::AtomProperty< long long int > & ( ::SireMol::AtomProperty< long long int >::*set_function_type )( ::SireMol::CGAtomIdx const &,long long int const & ) ;
            set_function_type set_function_value( &::SireMol::AtomProperty< long long int >::set );
            
            AtomIntProperty_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("cgatomidx"), bp::arg("value") )
                , bp::return_self< >() );
        
        }
        { //::SireMol::AtomProperty< long long int >::set
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef ::SireMol::AtomProperty< long long int > & ( ::SireMol::AtomProperty< long long int >::*set_function_type )( ::SireMol::CGIdx,::QVector< long long int > const & ) ;
            set_function_type set_function_value( &::SireMol::AtomProperty< long long int >::set );
            
            AtomIntProperty_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("cgidx"), bp::arg("values") )
                , bp::return_self< >() );
        
        }
        { //::SireMol::AtomProperty< long long int >::size
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< long long int >::*size_function_type )(  ) const;
            size_function_type size_function_value( &::SireMol::AtomProperty< long long int >::size );
            
            AtomIntProperty_exposer.def( 
                "size"
                , size_function_value );
        
        }
        { //::SireMol::AtomProperty< long long int >::toString
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef ::QString ( ::SireMol::AtomProperty< long long int >::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMol::AtomProperty< long long int >::toString );
            
            AtomIntProperty_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMol::AtomProperty< long long int >::toVariant
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef ::SireMol::AtomProperty< QVariant > ( ::SireMol::AtomProperty< long long int >::*toVariant_function_type )(  ) const;
            toVariant_function_type toVariant_function_value( &::SireMol::AtomProperty< long long int >::toVariant );
            
            AtomIntProperty_exposer.def( 
                "toVariant"
                , toVariant_function_value );
        
        }
        { //::SireMol::AtomProperty< long long int >::toVector
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef ::QVector< long long int > ( ::SireMol::AtomProperty< long long int >::*toVector_function_type )(  ) const;
            toVector_function_type toVector_function_value( &::SireMol::AtomProperty< long long int >::toVector );
            
            AtomIntProperty_exposer.def( 
                "toVector"
                , toVector_function_value );
        
        }
        { //::SireMol::AtomProperty< long long int >::toVector
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef ::QVector< long long int > ( ::SireMol::AtomProperty< long long int >::*toVector_function_type )( ::SireMol::AtomSelection const & ) const;
            toVector_function_type toVector_function_value( &::SireMol::AtomProperty< long long int >::toVector );
            
            AtomIntProperty_exposer.def( 
                "toVector"
                , toVector_function_value
                , ( bp::arg("selection") ) );
        
        }
        { //::SireMol::AtomProperty< long long int >::typeName
        
            typedef SireMol::AtomProperty< long long int > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::AtomProperty< long long int >::typeName );
            
            AtomIntProperty_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        AtomIntProperty_exposer.staticmethod( "fromVariant" );
        AtomIntProperty_exposer.staticmethod( "typeName" );
        AtomIntProperty_exposer.def( "__copy__", &__copy__);
        AtomIntProperty_exposer.def( "__deepcopy__", &__copy__);
        AtomIntProperty_exposer.def( "clone", &__copy__);
        AtomIntProperty_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::AtomProperty<long long> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomIntProperty_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::AtomProperty<long long> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomIntProperty_exposer.def( "__str__", &__str__< ::SireMol::AtomProperty<long long> > );
        AtomIntProperty_exposer.def( "__repr__", &__str__< ::SireMol::AtomProperty<long long> > );
    }

}
