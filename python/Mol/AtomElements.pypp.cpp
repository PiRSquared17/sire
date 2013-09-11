// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "AtomElements.pypp.hpp"

namespace bp = boost::python;

#include "atomelements.h"

#include "atomelements.h"

#include "SireMaths/vector.h"

SireMol::AtomProperty<SireMol::Element> __copy__(const SireMol::AtomProperty<SireMol::Element> &other){ return SireMol::AtomProperty<SireMol::Element>(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_AtomElements_class(){

    { //::SireMol::AtomProperty< SireMol::Element >
        typedef bp::class_< SireMol::AtomProperty< SireMol::Element >, bp::bases< SireMol::AtomProp, SireMol::MolViewProperty, SireBase::Property > > AtomElements_exposer_t;
        AtomElements_exposer_t AtomElements_exposer = AtomElements_exposer_t( "AtomElements", bp::init< >() );
        bp::scope AtomElements_scope( AtomElements_exposer );
        AtomElements_exposer.def( bp::init< SireMol::MoleculeInfoData const & >(( bp::arg("molinfo") )) );
        AtomElements_exposer.def( bp::init< SireMol::MoleculeInfoData const &, SireMol::Element const & >(( bp::arg("molinfo"), bp::arg("default_value") )) );
        AtomElements_exposer.def( bp::init< SireMol::Element const & >(( bp::arg("value") )) );
        AtomElements_exposer.def( bp::init< SireBase::PackedArray2D< SireMol::Element > const & >(( bp::arg("values") )) );
        AtomElements_exposer.def( bp::init< SireMol::AtomProperty< SireMol::Element > const & >(( bp::arg("other") )) );
        { //::SireMol::AtomProperty< SireMol::Element >::array
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef ::SireBase::PackedArray2D< SireMol::Element > const & ( ::SireMol::AtomProperty< SireMol::Element >::*array_function_type )(  ) const;
            array_function_type array_function_value( &::SireMol::AtomProperty< SireMol::Element >::array );
            
            AtomElements_exposer.def( 
                "array"
                , array_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::assertCanConvert
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< SireMol::Element >::*assertCanConvert_function_type )( ::QVariant const & ) const;
            assertCanConvert_function_type assertCanConvert_function_value( &::SireMol::AtomProperty< SireMol::Element >::assertCanConvert );
            
            AtomElements_exposer.def( 
                "assertCanConvert"
                , assertCanConvert_function_value
                , ( bp::arg("value") ) );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::assignFrom
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< SireMol::Element >::*assignFrom_function_type )( ::SireMol::AtomProperty< QVariant > const & ) ;
            assignFrom_function_type assignFrom_function_value( &::SireMol::AtomProperty< SireMol::Element >::assignFrom );
            
            AtomElements_exposer.def( 
                "assignFrom"
                , assignFrom_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::at
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef ::SireBase::detail::PackedArray2D_Array< SireMol::Element > const & ( ::SireMol::AtomProperty< SireMol::Element >::*at_function_type )( ::SireMol::CGIdx ) const;
            at_function_type at_function_value( &::SireMol::AtomProperty< SireMol::Element >::at );
            
            AtomElements_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("cgidx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::at
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef ::SireMol::Element const & ( ::SireMol::AtomProperty< SireMol::Element >::*at_function_type )( ::SireMol::CGAtomIdx const & ) const;
            at_function_type at_function_value( &::SireMol::AtomProperty< SireMol::Element >::at );
            
            AtomElements_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("cgatomidx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::canConvert
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< SireMol::Element >::*canConvert_function_type )( ::QVariant const & ) const;
            canConvert_function_type canConvert_function_value( &::SireMol::AtomProperty< SireMol::Element >::canConvert );
            
            AtomElements_exposer.def( 
                "canConvert"
                , canConvert_function_value
                , ( bp::arg("value") ) );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::copyFrom
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< SireMol::Element >::*copyFrom_function_type )( ::QVector< SireMol::Element > const & ) ;
            copyFrom_function_type copyFrom_function_value( &::SireMol::AtomProperty< SireMol::Element >::copyFrom );
            
            AtomElements_exposer.def( 
                "copyFrom"
                , copyFrom_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::copyFrom
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< SireMol::Element >::*copyFrom_function_type )( ::QVector< SireMol::Element > const &,::SireMol::AtomSelection const & ) ;
            copyFrom_function_type copyFrom_function_value( &::SireMol::AtomProperty< SireMol::Element >::copyFrom );
            
            AtomElements_exposer.def( 
                "copyFrom"
                , copyFrom_function_value
                , ( bp::arg("values"), bp::arg("selection") ) );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::count
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireMol::Element >::*count_function_type )(  ) const;
            count_function_type count_function_value( &::SireMol::AtomProperty< SireMol::Element >::count );
            
            AtomElements_exposer.def( 
                "count"
                , count_function_value );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::divide
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef ::SireBase::PropPtr< SireBase::Property > ( ::SireMol::AtomProperty< SireMol::Element >::*divide_function_type )( ::QVector< SireMol::AtomSelection > const & ) const;
            divide_function_type divide_function_value( &::SireMol::AtomProperty< SireMol::Element >::divide );
            
            AtomElements_exposer.def( 
                "divide"
                , divide_function_value
                , ( bp::arg("beads") ) );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::divideByResidue
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef ::SireBase::PropPtr< SireBase::Property > ( ::SireMol::AtomProperty< SireMol::Element >::*divideByResidue_function_type )( ::SireMol::MoleculeInfoData const & ) const;
            divideByResidue_function_type divideByResidue_function_value( &::SireMol::AtomProperty< SireMol::Element >::divideByResidue );
            
            AtomElements_exposer.def( 
                "divideByResidue"
                , divideByResidue_function_value
                , ( bp::arg("molinfo") ) );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::fromVariant
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef ::SireMol::AtomProperty< SireMol::Element > ( *fromVariant_function_type )( ::SireMol::AtomProperty< QVariant > const & );
            fromVariant_function_type fromVariant_function_value( &::SireMol::AtomProperty< SireMol::Element >::fromVariant );
            
            AtomElements_exposer.def( 
                "fromVariant"
                , fromVariant_function_value
                , ( bp::arg("variant") ) );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::get
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef ::SireBase::detail::PackedArray2D_Array< SireMol::Element > const & ( ::SireMol::AtomProperty< SireMol::Element >::*get_function_type )( ::SireMol::CGIdx ) const;
            get_function_type get_function_value( &::SireMol::AtomProperty< SireMol::Element >::get );
            
            AtomElements_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("cgidx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::get
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef ::SireMol::Element const & ( ::SireMol::AtomProperty< SireMol::Element >::*get_function_type )( ::SireMol::CGAtomIdx const & ) const;
            get_function_type get_function_value( &::SireMol::AtomProperty< SireMol::Element >::get );
            
            AtomElements_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("cgatomidx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::isCompatibleWith
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< SireMol::Element >::*isCompatibleWith_function_type )( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::AtomProperty< SireMol::Element >::isCompatibleWith );
            
            AtomElements_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") ) );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::isEmpty
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< SireMol::Element >::*isEmpty_function_type )(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::AtomProperty< SireMol::Element >::isEmpty );
            
            AtomElements_exposer.def( 
                "isEmpty"
                , isEmpty_function_value );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::matchToSelection
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef ::SireMol::AtomProperty< SireMol::Element > ( ::SireMol::AtomProperty< SireMol::Element >::*matchToSelection_function_type )( ::SireMol::AtomSelection const & ) const;
            matchToSelection_function_type matchToSelection_function_value( &::SireMol::AtomProperty< SireMol::Element >::matchToSelection );
            
            AtomElements_exposer.def( 
                "matchToSelection"
                , matchToSelection_function_value
                , ( bp::arg("selection") ) );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::merge
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef ::SireBase::PropPtr< SireBase::Property > ( ::SireMol::AtomProperty< SireMol::Element >::*merge_function_type )( ::SireMol::MoleculeInfoData const & ) const;
            merge_function_type merge_function_value( &::SireMol::AtomProperty< SireMol::Element >::merge );
            
            AtomElements_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("molinfo") ) );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::nAtoms
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireMol::Element >::*nAtoms_function_type )(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::AtomProperty< SireMol::Element >::nAtoms );
            
            AtomElements_exposer.def( 
                "nAtoms"
                , nAtoms_function_value );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::nAtoms
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireMol::Element >::*nAtoms_function_type )( ::SireMol::CGIdx ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::AtomProperty< SireMol::Element >::nAtoms );
            
            AtomElements_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , ( bp::arg("cgidx") ) );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::nCutGroups
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireMol::Element >::*nCutGroups_function_type )(  ) const;
            nCutGroups_function_type nCutGroups_function_value( &::SireMol::AtomProperty< SireMol::Element >::nCutGroups );
            
            AtomElements_exposer.def( 
                "nCutGroups"
                , nCutGroups_function_value );
        
        }
        AtomElements_exposer.def( bp::self != bp::self );
        { //::SireMol::AtomProperty< SireMol::Element >::operator=
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef ::SireMol::AtomProperty< SireMol::Element > & ( ::SireMol::AtomProperty< SireMol::Element >::*assign_function_type )( ::SireMol::AtomProperty< SireMol::Element > const & ) ;
            assign_function_type assign_function_value( &::SireMol::AtomProperty< SireMol::Element >::operator= );
            
            AtomElements_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        AtomElements_exposer.def( bp::self == bp::self );
        { //::SireMol::AtomProperty< SireMol::Element >::operator[]
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef ::SireBase::detail::PackedArray2D_Array< SireMol::Element > const & ( ::SireMol::AtomProperty< SireMol::Element >::*__getitem___function_type )( ::SireMol::CGIdx ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::AtomProperty< SireMol::Element >::operator[] );
            
            AtomElements_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("cgidx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::operator[]
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef ::SireMol::Element const & ( ::SireMol::AtomProperty< SireMol::Element >::*__getitem___function_type )( ::SireMol::CGAtomIdx const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::AtomProperty< SireMol::Element >::operator[] );
            
            AtomElements_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("cgatomidx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::set
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef ::SireMol::AtomProperty< SireMol::Element > & ( ::SireMol::AtomProperty< SireMol::Element >::*set_function_type )( ::SireMol::CGAtomIdx const &,::SireMol::Element const & ) ;
            set_function_type set_function_value( &::SireMol::AtomProperty< SireMol::Element >::set );
            
            AtomElements_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("cgatomidx"), bp::arg("value") )
                , bp::return_self< >() );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::set
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef ::SireMol::AtomProperty< SireMol::Element > & ( ::SireMol::AtomProperty< SireMol::Element >::*set_function_type )( ::SireMol::CGIdx,::QVector< SireMol::Element > const & ) ;
            set_function_type set_function_value( &::SireMol::AtomProperty< SireMol::Element >::set );
            
            AtomElements_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("cgidx"), bp::arg("values") )
                , bp::return_self< >() );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::size
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireMol::Element >::*size_function_type )(  ) const;
            size_function_type size_function_value( &::SireMol::AtomProperty< SireMol::Element >::size );
            
            AtomElements_exposer.def( 
                "size"
                , size_function_value );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::toString
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef ::QString ( ::SireMol::AtomProperty< SireMol::Element >::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMol::AtomProperty< SireMol::Element >::toString );
            
            AtomElements_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::toVariant
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef ::SireMol::AtomProperty< QVariant > ( ::SireMol::AtomProperty< SireMol::Element >::*toVariant_function_type )(  ) const;
            toVariant_function_type toVariant_function_value( &::SireMol::AtomProperty< SireMol::Element >::toVariant );
            
            AtomElements_exposer.def( 
                "toVariant"
                , toVariant_function_value );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::toVector
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef ::QVector< SireMol::Element > ( ::SireMol::AtomProperty< SireMol::Element >::*toVector_function_type )(  ) const;
            toVector_function_type toVector_function_value( &::SireMol::AtomProperty< SireMol::Element >::toVector );
            
            AtomElements_exposer.def( 
                "toVector"
                , toVector_function_value );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::toVector
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef ::QVector< SireMol::Element > ( ::SireMol::AtomProperty< SireMol::Element >::*toVector_function_type )( ::SireMol::AtomSelection const & ) const;
            toVector_function_type toVector_function_value( &::SireMol::AtomProperty< SireMol::Element >::toVector );
            
            AtomElements_exposer.def( 
                "toVector"
                , toVector_function_value
                , ( bp::arg("selection") ) );
        
        }
        { //::SireMol::AtomProperty< SireMol::Element >::typeName
        
            typedef SireMol::AtomProperty< SireMol::Element > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::AtomProperty< SireMol::Element >::typeName );
            
            AtomElements_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        AtomElements_exposer.staticmethod( "fromVariant" );
        AtomElements_exposer.staticmethod( "typeName" );
        AtomElements_exposer.def( "__copy__", &__copy__);
        AtomElements_exposer.def( "__deepcopy__", &__copy__);
        AtomElements_exposer.def( "clone", &__copy__);
        AtomElements_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::AtomProperty<SireMol::Element> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomElements_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::AtomProperty<SireMol::Element> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomElements_exposer.def( "__str__", &__str__< ::SireMol::AtomProperty<SireMol::Element> > );
        AtomElements_exposer.def( "__repr__", &__str__< ::SireMol::AtomProperty<SireMol::Element> > );
    }

}
