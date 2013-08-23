// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "CutGroup.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atom.h"

#include "cgatomidx.h"

#include "cgeditor.h"

#include "cutgroup.h"

#include "evaluator.h"

#include "groupatomids.h"

#include "molecule.h"

#include "mover.hpp"

#include "mover_metaid.h"

#include "selector.hpp"

#include <QDebug>

#include "cutgroup.h"

#include "cgproperty.hpp"

const QString& get_Metadata_SireMol_CGStringProperty_function1(const SireMol::CutGroup &atom,
                                   const QString &metakey){ return atom.metadata< QString >(metakey); }

const QString& get_Metadata_SireMol_CGStringProperty_function2(const SireMol::CutGroup &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< QString >(key, metakey); }

const qint64& get_Metadata_SireMol_CGIntProperty_function1(const SireMol::CutGroup &atom,
                                   const QString &metakey){ return atom.metadata< qint64 >(metakey); }

const qint64& get_Metadata_SireMol_CGIntProperty_function2(const SireMol::CutGroup &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< qint64 >(key, metakey); }

const double& get_Metadata_SireMol_CGFloatProperty_function1(const SireMol::CutGroup &atom,
                                   const QString &metakey){ return atom.metadata< double >(metakey); }

const double& get_Metadata_SireMol_CGFloatProperty_function2(const SireMol::CutGroup &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< double >(key, metakey); }

const QVariant& get_Metadata_SireMol_CGVariantProperty_function1(const SireMol::CutGroup &atom,
                                   const QString &metakey){ return atom.metadata< QVariant >(metakey); }

const QVariant& get_Metadata_SireMol_CGVariantProperty_function2(const SireMol::CutGroup &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< QVariant >(key, metakey); }

SireMol::CutGroup __copy__(const SireMol::CutGroup &other){ return SireMol::CutGroup(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_CutGroup_class(){

    { //::SireMol::CutGroup
        typedef bp::class_< SireMol::CutGroup, bp::bases< SireMol::MoleculeView, SireBase::Property > > CutGroup_exposer_t;
        CutGroup_exposer_t CutGroup_exposer = CutGroup_exposer_t( "CutGroup", bp::init< >() );
        bp::scope CutGroup_scope( CutGroup_exposer );
        CutGroup_exposer.def( bp::init< SireMol::MoleculeData const &, SireMol::CGID const & >(( bp::arg("moldata"), bp::arg("cgid") )) );
        CutGroup_exposer.def( bp::init< SireMol::CutGroup const & >(( bp::arg("other") )) );
        { //::SireMol::CutGroup::assertContainsMetadata
        
            typedef void ( ::SireMol::CutGroup::*assertContainsMetadata_function_type )( ::SireBase::PropertyName const & ) const;
            assertContainsMetadata_function_type assertContainsMetadata_function_value( &::SireMol::CutGroup::assertContainsMetadata );
            
            CutGroup_exposer.def( 
                "assertContainsMetadata"
                , assertContainsMetadata_function_value
                , ( bp::arg("metakey") ) );
        
        }
        { //::SireMol::CutGroup::assertContainsMetadata
        
            typedef void ( ::SireMol::CutGroup::*assertContainsMetadata_function_type )( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            assertContainsMetadata_function_type assertContainsMetadata_function_value( &::SireMol::CutGroup::assertContainsMetadata );
            
            CutGroup_exposer.def( 
                "assertContainsMetadata"
                , assertContainsMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") ) );
        
        }
        { //::SireMol::CutGroup::assertContainsProperty
        
            typedef void ( ::SireMol::CutGroup::*assertContainsProperty_function_type )( ::SireBase::PropertyName const & ) const;
            assertContainsProperty_function_type assertContainsProperty_function_value( &::SireMol::CutGroup::assertContainsProperty );
            
            CutGroup_exposer.def( 
                "assertContainsProperty"
                , assertContainsProperty_function_value
                , ( bp::arg("key") ) );
        
        }
        { //::SireMol::CutGroup::atomIdxs
        
            typedef ::QList< SireMol::AtomIdx > const & ( ::SireMol::CutGroup::*atomIdxs_function_type )(  ) const;
            atomIdxs_function_type atomIdxs_function_value( &::SireMol::CutGroup::atomIdxs );
            
            CutGroup_exposer.def( 
                "atomIdxs"
                , atomIdxs_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::CutGroup::contains
        
            typedef bool ( ::SireMol::CutGroup::*contains_function_type )( ::SireMol::AtomID const & ) const;
            contains_function_type contains_function_value( &::SireMol::CutGroup::contains );
            
            CutGroup_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("atomid") ) );
        
        }
        { //::SireMol::CutGroup::contains
        
            typedef bool ( ::SireMol::CutGroup::*contains_function_type )( ::SireMol::AtomIdx ) const;
            contains_function_type contains_function_value( &::SireMol::CutGroup::contains );
            
            CutGroup_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("atomidx") ) );
        
        }
        { //::SireMol::CutGroup::edit
        
            typedef ::SireMol::CGEditor ( ::SireMol::CutGroup::*edit_function_type )(  ) const;
            edit_function_type edit_function_value( &::SireMol::CutGroup::edit );
            
            CutGroup_exposer.def( 
                "edit"
                , edit_function_value );
        
        }
        { //::SireMol::CutGroup::evaluate
        
            typedef ::SireMol::Evaluator ( ::SireMol::CutGroup::*evaluate_function_type )(  ) const;
            evaluate_function_type evaluate_function_value( &::SireMol::CutGroup::evaluate );
            
            CutGroup_exposer.def( 
                "evaluate"
                , evaluate_function_value );
        
        }
        { //::SireMol::CutGroup::hasMetadata
        
            typedef bool ( ::SireMol::CutGroup::*hasMetadata_function_type )( ::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::CutGroup::hasMetadata );
            
            CutGroup_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("metakey") ) );
        
        }
        { //::SireMol::CutGroup::hasMetadata
        
            typedef bool ( ::SireMol::CutGroup::*hasMetadata_function_type )( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::CutGroup::hasMetadata );
            
            CutGroup_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") ) );
        
        }
        { //::SireMol::CutGroup::hasProperty
        
            typedef bool ( ::SireMol::CutGroup::*hasProperty_function_type )( ::SireBase::PropertyName const & ) const;
            hasProperty_function_type hasProperty_function_value( &::SireMol::CutGroup::hasProperty );
            
            CutGroup_exposer.def( 
                "hasProperty"
                , hasProperty_function_value
                , ( bp::arg("key") ) );
        
        }
        { //::SireMol::CutGroup::index
        
            typedef ::SireMol::CGIdx ( ::SireMol::CutGroup::*index_function_type )(  ) const;
            index_function_type index_function_value( &::SireMol::CutGroup::index );
            
            CutGroup_exposer.def( 
                "index"
                , index_function_value );
        
        }
        { //::SireMol::CutGroup::intersects
        
            typedef bool ( ::SireMol::CutGroup::*intersects_function_type )( ::SireMol::AtomID const & ) const;
            intersects_function_type intersects_function_value( &::SireMol::CutGroup::intersects );
            
            CutGroup_exposer.def( 
                "intersects"
                , intersects_function_value
                , ( bp::arg("atomid") ) );
        
        }
        { //::SireMol::CutGroup::isEmpty
        
            typedef bool ( ::SireMol::CutGroup::*isEmpty_function_type )(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::CutGroup::isEmpty );
            
            CutGroup_exposer.def( 
                "isEmpty"
                , isEmpty_function_value );
        
        }
        { //::SireMol::CutGroup::metadataKeys
        
            typedef ::QStringList ( ::SireMol::CutGroup::*metadataKeys_function_type )(  ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::CutGroup::metadataKeys );
            
            CutGroup_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value );
        
        }
        { //::SireMol::CutGroup::metadataKeys
        
            typedef ::QStringList ( ::SireMol::CutGroup::*metadataKeys_function_type )( ::SireBase::PropertyName const & ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::CutGroup::metadataKeys );
            
            CutGroup_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , ( bp::arg("key") ) );
        
        }
        { //::SireMol::CutGroup::move
        
            typedef ::SireMol::Mover< SireMol::CutGroup > ( ::SireMol::CutGroup::*move_function_type )(  ) const;
            move_function_type move_function_value( &::SireMol::CutGroup::move );
            
            CutGroup_exposer.def( 
                "move"
                , move_function_value );
        
        }
        { //::SireMol::CutGroup::nAtoms
        
            typedef int ( ::SireMol::CutGroup::*nAtoms_function_type )(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::CutGroup::nAtoms );
            
            CutGroup_exposer.def( 
                "nAtoms"
                , nAtoms_function_value );
        
        }
        { //::SireMol::CutGroup::name
        
            typedef ::SireMol::CGName const & ( ::SireMol::CutGroup::*name_function_type )(  ) const;
            name_function_type name_function_value( &::SireMol::CutGroup::name );
            
            CutGroup_exposer.def( 
                "name"
                , name_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        CutGroup_exposer.def( bp::self != bp::self );
        { //::SireMol::CutGroup::operator=
        
            typedef ::SireMol::CutGroup & ( ::SireMol::CutGroup::*assign_function_type )( ::SireMol::CutGroup const & ) ;
            assign_function_type assign_function_value( &::SireMol::CutGroup::operator= );
            
            CutGroup_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        CutGroup_exposer.def( bp::self == bp::self );
        { //::SireMol::CutGroup::propertyKeys
        
            typedef ::QStringList ( ::SireMol::CutGroup::*propertyKeys_function_type )(  ) const;
            propertyKeys_function_type propertyKeys_function_value( &::SireMol::CutGroup::propertyKeys );
            
            CutGroup_exposer.def( 
                "propertyKeys"
                , propertyKeys_function_value );
        
        }
        { //::SireMol::CutGroup::selectedAll
        
            typedef bool ( ::SireMol::CutGroup::*selectedAll_function_type )(  ) const;
            selectedAll_function_type selectedAll_function_value( &::SireMol::CutGroup::selectedAll );
            
            CutGroup_exposer.def( 
                "selectedAll"
                , selectedAll_function_value );
        
        }
        { //::SireMol::CutGroup::selection
        
            typedef ::SireMol::AtomSelection ( ::SireMol::CutGroup::*selection_function_type )(  ) const;
            selection_function_type selection_function_value( &::SireMol::CutGroup::selection );
            
            CutGroup_exposer.def( 
                "selection"
                , selection_function_value );
        
        }
        { //::SireMol::CutGroup::selector
        
            typedef ::SireMol::Selector< SireMol::CutGroup > ( ::SireMol::CutGroup::*selector_function_type )(  ) const;
            selector_function_type selector_function_value( &::SireMol::CutGroup::selector );
            
            CutGroup_exposer.def( 
                "selector"
                , selector_function_value );
        
        }
        { //::SireMol::CutGroup::toString
        
            typedef ::QString ( ::SireMol::CutGroup::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMol::CutGroup::toString );
            
            CutGroup_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMol::CutGroup::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::CutGroup::typeName );
            
            CutGroup_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMol::CutGroup::update
        
            typedef void ( ::SireMol::CutGroup::*update_function_type )( ::SireMol::MoleculeData const & ) ;
            update_function_type update_function_value( &::SireMol::CutGroup::update );
            
            CutGroup_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("moldata") ) );
        
        }
        CutGroup_exposer.staticmethod( "typeName" );
        CutGroup_exposer.def( "_get_property_SireMol_CGStringProperty", &SireMol::CutGroup::property< QString >, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_metadata_SireMol_CGStringProperty", get_Metadata_SireMol_CGStringProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_metadata_SireMol_CGStringProperty", &get_Metadata_SireMol_CGStringProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_property_SireMol_CGIntProperty", &SireMol::CutGroup::property< qint64 >, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_metadata_SireMol_CGIntProperty", get_Metadata_SireMol_CGIntProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_metadata_SireMol_CGIntProperty", &get_Metadata_SireMol_CGIntProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_property_SireMol_CGFloatProperty", &SireMol::CutGroup::property< double >, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_metadata_SireMol_CGFloatProperty", get_Metadata_SireMol_CGFloatProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_metadata_SireMol_CGFloatProperty", &get_Metadata_SireMol_CGFloatProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_property_SireMol_CGVariantProperty", &SireMol::CutGroup::property< QVariant >, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_metadata_SireMol_CGVariantProperty", get_Metadata_SireMol_CGVariantProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_metadata_SireMol_CGVariantProperty", &get_Metadata_SireMol_CGVariantProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "__copy__", &__copy__);
        CutGroup_exposer.def( "__deepcopy__", &__copy__);
        CutGroup_exposer.def( "clone", &__copy__);
        CutGroup_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::CutGroup >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CutGroup_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::CutGroup >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CutGroup_exposer.def( "__str__", &__str__< ::SireMol::CutGroup > );
        CutGroup_exposer.def( "__repr__", &__str__< ::SireMol::CutGroup > );
    }

}