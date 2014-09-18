// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Residue.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireError/errors.h"

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atom.h"

#include "atomidcombos.h"

#include "chain.h"

#include "evaluator.h"

#include "molecule.h"

#include "mover.hpp"

#include "mover_metaid.h"

#include "reseditor.h"

#include "residue.h"

#include "selector.hpp"

#include "residue.h"

#include "resproperty.hpp"

const QString& get_Metadata_SireMol_ResStringProperty_function1(const SireMol::Residue &atom,
                                   const QString &metakey){ return atom.metadata< QString >(metakey); }

const QString& get_Metadata_SireMol_ResStringProperty_function2(const SireMol::Residue &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< QString >(key, metakey); }

const qint64& get_Metadata_SireMol_ResIntProperty_function1(const SireMol::Residue &atom,
                                   const QString &metakey){ return atom.metadata< qint64 >(metakey); }

const qint64& get_Metadata_SireMol_ResIntProperty_function2(const SireMol::Residue &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< qint64 >(key, metakey); }

const double& get_Metadata_SireMol_ResFloatProperty_function1(const SireMol::Residue &atom,
                                   const QString &metakey){ return atom.metadata< double >(metakey); }

const double& get_Metadata_SireMol_ResFloatProperty_function2(const SireMol::Residue &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< double >(key, metakey); }

const QVariant& get_Metadata_SireMol_ResVariantProperty_function1(const SireMol::Residue &atom,
                                   const QString &metakey){ return atom.metadata< QVariant >(metakey); }

const QVariant& get_Metadata_SireMol_ResVariantProperty_function2(const SireMol::Residue &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< QVariant >(key, metakey); }

SireMol::Residue __copy__(const SireMol::Residue &other){ return SireMol::Residue(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_Residue_class(){

    { //::SireMol::Residue
        typedef bp::class_< SireMol::Residue, bp::bases< SireMol::MoleculeView, SireBase::Property > > Residue_exposer_t;
        Residue_exposer_t Residue_exposer = Residue_exposer_t( "Residue", bp::init< >() );
        bp::scope Residue_scope( Residue_exposer );
        Residue_exposer.def( bp::init< SireMol::MoleculeData const &, SireMol::ResID const & >(( bp::arg("moldata"), bp::arg("resid") )) );
        Residue_exposer.def( bp::init< SireMol::Residue const & >(( bp::arg("other") )) );
        { //::SireMol::Residue::assertContainsMetadata
        
            typedef void ( ::SireMol::Residue::*assertContainsMetadata_function_type )( ::SireBase::PropertyName const & ) const;
            assertContainsMetadata_function_type assertContainsMetadata_function_value( &::SireMol::Residue::assertContainsMetadata );
            
            Residue_exposer.def( 
                "assertContainsMetadata"
                , assertContainsMetadata_function_value
                , ( bp::arg("metakey") ) );
        
        }
        { //::SireMol::Residue::assertContainsMetadata
        
            typedef void ( ::SireMol::Residue::*assertContainsMetadata_function_type )( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            assertContainsMetadata_function_type assertContainsMetadata_function_value( &::SireMol::Residue::assertContainsMetadata );
            
            Residue_exposer.def( 
                "assertContainsMetadata"
                , assertContainsMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") ) );
        
        }
        { //::SireMol::Residue::assertContainsProperty
        
            typedef void ( ::SireMol::Residue::*assertContainsProperty_function_type )( ::SireBase::PropertyName const & ) const;
            assertContainsProperty_function_type assertContainsProperty_function_value( &::SireMol::Residue::assertContainsProperty );
            
            Residue_exposer.def( 
                "assertContainsProperty"
                , assertContainsProperty_function_value
                , ( bp::arg("key") ) );
        
        }
        { //::SireMol::Residue::atomIdxs
        
            typedef ::QList< SireMol::AtomIdx > const & ( ::SireMol::Residue::*atomIdxs_function_type )(  ) const;
            atomIdxs_function_type atomIdxs_function_value( &::SireMol::Residue::atomIdxs );
            
            Residue_exposer.def( 
                "atomIdxs"
                , atomIdxs_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::Residue::contains
        
            typedef bool ( ::SireMol::Residue::*contains_function_type )( ::SireMol::AtomIdx ) const;
            contains_function_type contains_function_value( &::SireMol::Residue::contains );
            
            Residue_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("atomidx") ) );
        
        }
        { //::SireMol::Residue::contains
        
            typedef bool ( ::SireMol::Residue::*contains_function_type )( ::SireMol::AtomID const & ) const;
            contains_function_type contains_function_value( &::SireMol::Residue::contains );
            
            Residue_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("atomid") ) );
        
        }
        { //::SireMol::Residue::edit
        
            typedef ::SireMol::ResEditor ( ::SireMol::Residue::*edit_function_type )(  ) const;
            edit_function_type edit_function_value( &::SireMol::Residue::edit );
            
            Residue_exposer.def( 
                "edit"
                , edit_function_value );
        
        }
        { //::SireMol::Residue::evaluate
        
            typedef ::SireMol::Evaluator ( ::SireMol::Residue::*evaluate_function_type )(  ) const;
            evaluate_function_type evaluate_function_value( &::SireMol::Residue::evaluate );
            
            Residue_exposer.def( 
                "evaluate"
                , evaluate_function_value );
        
        }
        { //::SireMol::Residue::hasMetadata
        
            typedef bool ( ::SireMol::Residue::*hasMetadata_function_type )( ::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::Residue::hasMetadata );
            
            Residue_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("metakey") ) );
        
        }
        { //::SireMol::Residue::hasMetadata
        
            typedef bool ( ::SireMol::Residue::*hasMetadata_function_type )( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::Residue::hasMetadata );
            
            Residue_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") ) );
        
        }
        { //::SireMol::Residue::hasProperty
        
            typedef bool ( ::SireMol::Residue::*hasProperty_function_type )( ::SireBase::PropertyName const & ) const;
            hasProperty_function_type hasProperty_function_value( &::SireMol::Residue::hasProperty );
            
            Residue_exposer.def( 
                "hasProperty"
                , hasProperty_function_value
                , ( bp::arg("key") ) );
        
        }
        { //::SireMol::Residue::index
        
            typedef ::SireMol::ResIdx ( ::SireMol::Residue::*index_function_type )(  ) const;
            index_function_type index_function_value( &::SireMol::Residue::index );
            
            Residue_exposer.def( 
                "index"
                , index_function_value );
        
        }
        { //::SireMol::Residue::intersects
        
            typedef bool ( ::SireMol::Residue::*intersects_function_type )( ::SireMol::AtomID const & ) const;
            intersects_function_type intersects_function_value( &::SireMol::Residue::intersects );
            
            Residue_exposer.def( 
                "intersects"
                , intersects_function_value
                , ( bp::arg("atomid") ) );
        
        }
        { //::SireMol::Residue::isEmpty
        
            typedef bool ( ::SireMol::Residue::*isEmpty_function_type )(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::Residue::isEmpty );
            
            Residue_exposer.def( 
                "isEmpty"
                , isEmpty_function_value );
        
        }
        { //::SireMol::Residue::isWithinChain
        
            typedef bool ( ::SireMol::Residue::*isWithinChain_function_type )(  ) const;
            isWithinChain_function_type isWithinChain_function_value( &::SireMol::Residue::isWithinChain );
            
            Residue_exposer.def( 
                "isWithinChain"
                , isWithinChain_function_value );
        
        }
        { //::SireMol::Residue::metadataKeys
        
            typedef ::QStringList ( ::SireMol::Residue::*metadataKeys_function_type )(  ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::Residue::metadataKeys );
            
            Residue_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value );
        
        }
        { //::SireMol::Residue::metadataKeys
        
            typedef ::QStringList ( ::SireMol::Residue::*metadataKeys_function_type )( ::SireBase::PropertyName const & ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::Residue::metadataKeys );
            
            Residue_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , ( bp::arg("key") ) );
        
        }
        { //::SireMol::Residue::move
        
            typedef ::SireMol::Mover< SireMol::Residue > ( ::SireMol::Residue::*move_function_type )(  ) const;
            move_function_type move_function_value( &::SireMol::Residue::move );
            
            Residue_exposer.def( 
                "move"
                , move_function_value );
        
        }
        { //::SireMol::Residue::nAtoms
        
            typedef int ( ::SireMol::Residue::*nAtoms_function_type )(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::Residue::nAtoms );
            
            Residue_exposer.def( 
                "nAtoms"
                , nAtoms_function_value );
        
        }
        { //::SireMol::Residue::name
        
            typedef ::SireMol::ResName ( ::SireMol::Residue::*name_function_type )(  ) const;
            name_function_type name_function_value( &::SireMol::Residue::name );
            
            Residue_exposer.def( 
                "name"
                , name_function_value );
        
        }
        { //::SireMol::Residue::number
        
            typedef ::SireMol::ResNum ( ::SireMol::Residue::*number_function_type )(  ) const;
            number_function_type number_function_value( &::SireMol::Residue::number );
            
            Residue_exposer.def( 
                "number"
                , number_function_value );
        
        }
        Residue_exposer.def( bp::self != bp::self );
        { //::SireMol::Residue::operator=
        
            typedef ::SireMol::Residue & ( ::SireMol::Residue::*assign_function_type )( ::SireMol::Residue const & ) ;
            assign_function_type assign_function_value( &::SireMol::Residue::operator= );
            
            Residue_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        Residue_exposer.def( bp::self == bp::self );
        { //::SireMol::Residue::propertyKeys
        
            typedef ::QStringList ( ::SireMol::Residue::*propertyKeys_function_type )(  ) const;
            propertyKeys_function_type propertyKeys_function_value( &::SireMol::Residue::propertyKeys );
            
            Residue_exposer.def( 
                "propertyKeys"
                , propertyKeys_function_value );
        
        }
        { //::SireMol::Residue::selectedAll
        
            typedef bool ( ::SireMol::Residue::*selectedAll_function_type )(  ) const;
            selectedAll_function_type selectedAll_function_value( &::SireMol::Residue::selectedAll );
            
            Residue_exposer.def( 
                "selectedAll"
                , selectedAll_function_value );
        
        }
        { //::SireMol::Residue::selection
        
            typedef ::SireMol::AtomSelection ( ::SireMol::Residue::*selection_function_type )(  ) const;
            selection_function_type selection_function_value( &::SireMol::Residue::selection );
            
            Residue_exposer.def( 
                "selection"
                , selection_function_value );
        
        }
        { //::SireMol::Residue::selector
        
            typedef ::SireMol::Selector< SireMol::Residue > ( ::SireMol::Residue::*selector_function_type )(  ) const;
            selector_function_type selector_function_value( &::SireMol::Residue::selector );
            
            Residue_exposer.def( 
                "selector"
                , selector_function_value );
        
        }
        { //::SireMol::Residue::toString
        
            typedef ::QString ( ::SireMol::Residue::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMol::Residue::toString );
            
            Residue_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMol::Residue::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::Residue::typeName );
            
            Residue_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMol::Residue::update
        
            typedef void ( ::SireMol::Residue::*update_function_type )( ::SireMol::MoleculeData const & ) ;
            update_function_type update_function_value( &::SireMol::Residue::update );
            
            Residue_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("moldata") ) );
        
        }
        Residue_exposer.staticmethod( "typeName" );
        Residue_exposer.def( "_get_property_SireMol_ResStringProperty", &SireMol::Residue::property< QString >, bp::return_value_policy<bp::copy_const_reference>());
        Residue_exposer.def( "_get_metadata_SireMol_ResStringProperty", get_Metadata_SireMol_ResStringProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        Residue_exposer.def( "_get_metadata_SireMol_ResStringProperty", &get_Metadata_SireMol_ResStringProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        Residue_exposer.def( "_get_property_SireMol_ResIntProperty", &SireMol::Residue::property< qint64 >, bp::return_value_policy<bp::copy_const_reference>());
        Residue_exposer.def( "_get_metadata_SireMol_ResIntProperty", get_Metadata_SireMol_ResIntProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        Residue_exposer.def( "_get_metadata_SireMol_ResIntProperty", &get_Metadata_SireMol_ResIntProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        Residue_exposer.def( "_get_property_SireMol_ResFloatProperty", &SireMol::Residue::property< double >, bp::return_value_policy<bp::copy_const_reference>());
        Residue_exposer.def( "_get_metadata_SireMol_ResFloatProperty", get_Metadata_SireMol_ResFloatProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        Residue_exposer.def( "_get_metadata_SireMol_ResFloatProperty", &get_Metadata_SireMol_ResFloatProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        Residue_exposer.def( "_get_property_SireMol_ResVariantProperty", &SireMol::Residue::property< QVariant >, bp::return_value_policy<bp::copy_const_reference>());
        Residue_exposer.def( "_get_metadata_SireMol_ResVariantProperty", get_Metadata_SireMol_ResVariantProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        Residue_exposer.def( "_get_metadata_SireMol_ResVariantProperty", &get_Metadata_SireMol_ResVariantProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        Residue_exposer.def( "__copy__", &__copy__);
        Residue_exposer.def( "__deepcopy__", &__copy__);
        Residue_exposer.def( "clone", &__copy__);
        Residue_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::Residue >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Residue_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::Residue >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Residue_exposer.def( "__str__", &__str__< ::SireMol::Residue > );
        Residue_exposer.def( "__repr__", &__str__< ::SireMol::Residue > );
    }

}
