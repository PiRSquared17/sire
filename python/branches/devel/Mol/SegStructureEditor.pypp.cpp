// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "SegStructureEditor.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "atom.h"

#include "atomeditor.h"

#include "chain.h"

#include "chaineditor.h"

#include "groupatomids.h"

#include "molecule.h"

#include "moleditor.h"

#include "mover.hpp"

#include "reseditor.h"

#include "residue.h"

#include "segeditor.h"

#include "segment.h"

#include "selector.hpp"

#include "segeditor.h"

SireMol::SegStructureEditor __copy__(const SireMol::SegStructureEditor &other){ return SireMol::SegStructureEditor(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_SegStructureEditor_class(){

    { //::SireMol::SegStructureEditor
        typedef bp::class_< SireMol::SegStructureEditor > SegStructureEditor_exposer_t;
        SegStructureEditor_exposer_t SegStructureEditor_exposer = SegStructureEditor_exposer_t( "SegStructureEditor", bp::init< >() );
        bp::scope SegStructureEditor_scope( SegStructureEditor_exposer );
        SegStructureEditor_exposer.def( bp::init< SireMol::Segment const & >(( bp::arg("residue") )) );
        SegStructureEditor_exposer.def( bp::init< SireMol::StructureEditor const &, SireMol::SegIdx >(( bp::arg("data"), bp::arg("residx") )) );
        SegStructureEditor_exposer.def( bp::init< SireMol::SegStructureEditor const & >(( bp::arg("other") )) );
        { //::SireMol::SegStructureEditor::add
        
            typedef ::SireMol::AtomStructureEditor ( ::SireMol::SegStructureEditor::*add_function_type )( ::SireMol::AtomName const & ) ;
            add_function_type add_function_value( &::SireMol::SegStructureEditor::add );
            
            SegStructureEditor_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("atomname") ) );
        
        }
        { //::SireMol::SegStructureEditor::add
        
            typedef ::SireMol::AtomStructureEditor ( ::SireMol::SegStructureEditor::*add_function_type )( ::SireMol::AtomNum ) ;
            add_function_type add_function_value( &::SireMol::SegStructureEditor::add );
            
            SegStructureEditor_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("atomnum") ) );
        
        }
        { //::SireMol::SegStructureEditor::atom
        
            typedef ::SireMol::AtomStructureEditor ( ::SireMol::SegStructureEditor::*atom_function_type )( int ) ;
            atom_function_type atom_function_value( &::SireMol::SegStructureEditor::atom );
            
            SegStructureEditor_exposer.def( 
                "atom"
                , atom_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireMol::SegStructureEditor::atom
        
            typedef ::SireMol::AtomStructureEditor ( ::SireMol::SegStructureEditor::*atom_function_type )( ::SireMol::AtomID const & ) ;
            atom_function_type atom_function_value( &::SireMol::SegStructureEditor::atom );
            
            SegStructureEditor_exposer.def( 
                "atom"
                , atom_function_value
                , ( bp::arg("atomid") ) );
        
        }
        { //::SireMol::SegStructureEditor::commit
        
            typedef ::SireMol::Segment ( ::SireMol::SegStructureEditor::*commit_function_type )(  ) const;
            commit_function_type commit_function_value( &::SireMol::SegStructureEditor::commit );
            
            SegStructureEditor_exposer.def( 
                "commit"
                , commit_function_value );
        
        }
        { //::SireMol::SegStructureEditor::index
        
            typedef ::SireMol::SegIdx ( ::SireMol::SegStructureEditor::*index_function_type )(  ) const;
            index_function_type index_function_value( &::SireMol::SegStructureEditor::index );
            
            SegStructureEditor_exposer.def( 
                "index"
                , index_function_value );
        
        }
        { //::SireMol::SegStructureEditor::molecule
        
            typedef ::SireMol::MolStructureEditor ( ::SireMol::SegStructureEditor::*molecule_function_type )(  ) ;
            molecule_function_type molecule_function_value( &::SireMol::SegStructureEditor::molecule );
            
            SegStructureEditor_exposer.def( 
                "molecule"
                , molecule_function_value );
        
        }
        { //::SireMol::SegStructureEditor::nAtoms
        
            typedef int ( ::SireMol::SegStructureEditor::*nAtoms_function_type )(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::SegStructureEditor::nAtoms );
            
            SegStructureEditor_exposer.def( 
                "nAtoms"
                , nAtoms_function_value );
        
        }
        { //::SireMol::SegStructureEditor::name
        
            typedef ::SireMol::SegName const & ( ::SireMol::SegStructureEditor::*name_function_type )(  ) const;
            name_function_type name_function_value( &::SireMol::SegStructureEditor::name );
            
            SegStructureEditor_exposer.def( 
                "name"
                , name_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMol::SegStructureEditor::operator=
        
            typedef ::SireMol::SegStructureEditor & ( ::SireMol::SegStructureEditor::*assign_function_type )( ::SireMol::Segment const & ) ;
            assign_function_type assign_function_value( &::SireMol::SegStructureEditor::operator= );
            
            SegStructureEditor_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("residue") )
                , bp::return_self< >() );
        
        }
        { //::SireMol::SegStructureEditor::operator=
        
            typedef ::SireMol::SegStructureEditor & ( ::SireMol::SegStructureEditor::*assign_function_type )( ::SireMol::SegStructureEditor const & ) ;
            assign_function_type assign_function_value( &::SireMol::SegStructureEditor::operator= );
            
            SegStructureEditor_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        { //::SireMol::SegStructureEditor::reindex
        
            typedef ::SireMol::SegStructureEditor & ( ::SireMol::SegStructureEditor::*reindex_function_type )( ::SireMol::SegIdx ) ;
            reindex_function_type reindex_function_value( &::SireMol::SegStructureEditor::reindex );
            
            SegStructureEditor_exposer.def( 
                "reindex"
                , reindex_function_value
                , ( bp::arg("index") )
                , bp::return_self< >() );
        
        }
        { //::SireMol::SegStructureEditor::remove
        
            typedef ::SireMol::MolStructureEditor ( ::SireMol::SegStructureEditor::*remove_function_type )(  ) ;
            remove_function_type remove_function_value( &::SireMol::SegStructureEditor::remove );
            
            SegStructureEditor_exposer.def( 
                "remove"
                , remove_function_value
                , bp::return_self< >() );
        
        }
        { //::SireMol::SegStructureEditor::remove
        
            typedef ::SireMol::SegStructureEditor & ( ::SireMol::SegStructureEditor::*remove_function_type )( ::SireMol::AtomID const & ) ;
            remove_function_type remove_function_value( &::SireMol::SegStructureEditor::remove );
            
            SegStructureEditor_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("atomid") )
                , bp::return_self< >() );
        
        }
        { //::SireMol::SegStructureEditor::remove
        
            typedef ::SireMol::SegStructureEditor & ( ::SireMol::SegStructureEditor::*remove_function_type )( int ) ;
            remove_function_type remove_function_value( &::SireMol::SegStructureEditor::remove );
            
            SegStructureEditor_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("i") )
                , bp::return_self< >() );
        
        }
        { //::SireMol::SegStructureEditor::rename
        
            typedef ::SireMol::SegStructureEditor & ( ::SireMol::SegStructureEditor::*rename_function_type )( ::SireMol::SegName const & ) ;
            rename_function_type rename_function_value( &::SireMol::SegStructureEditor::rename );
            
            SegStructureEditor_exposer.def( 
                "rename"
                , rename_function_value
                , ( bp::arg("name") )
                , bp::return_self< >() );
        
        }
        { //::SireMol::SegStructureEditor::select
        
            typedef ::SireMol::AtomStructureEditor ( ::SireMol::SegStructureEditor::*select_function_type )( int ) ;
            select_function_type select_function_value( &::SireMol::SegStructureEditor::select );
            
            SegStructureEditor_exposer.def( 
                "select"
                , select_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireMol::SegStructureEditor::select
        
            typedef ::SireMol::AtomStructureEditor ( ::SireMol::SegStructureEditor::*select_function_type )( ::SireMol::AtomID const & ) ;
            select_function_type select_function_value( &::SireMol::SegStructureEditor::select );
            
            SegStructureEditor_exposer.def( 
                "select"
                , select_function_value
                , ( bp::arg("atomid") ) );
        
        }
        { //::SireMol::SegStructureEditor::selectedAll
        
            typedef bool ( ::SireMol::SegStructureEditor::*selectedAll_function_type )(  ) const;
            selectedAll_function_type selectedAll_function_value( &::SireMol::SegStructureEditor::selectedAll );
            
            SegStructureEditor_exposer.def( 
                "selectedAll"
                , selectedAll_function_value );
        
        }
        { //::SireMol::SegStructureEditor::toString
        
            typedef ::QString ( ::SireMol::SegStructureEditor::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMol::SegStructureEditor::toString );
            
            SegStructureEditor_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMol::SegStructureEditor::transfer
        
            typedef ::SireMol::SegStructureEditor & ( ::SireMol::SegStructureEditor::*transfer_function_type )( ::SireMol::AtomID const &,::SireMol::SegID const & ) ;
            transfer_function_type transfer_function_value( &::SireMol::SegStructureEditor::transfer );
            
            SegStructureEditor_exposer.def( 
                "transfer"
                , transfer_function_value
                , ( bp::arg("atomid"), bp::arg("segid") )
                , bp::return_self< >() );
        
        }
        { //::SireMol::SegStructureEditor::transfer
        
            typedef ::SireMol::SegStructureEditor & ( ::SireMol::SegStructureEditor::*transfer_function_type )( int,::SireMol::SegID const & ) ;
            transfer_function_type transfer_function_value( &::SireMol::SegStructureEditor::transfer );
            
            SegStructureEditor_exposer.def( 
                "transfer"
                , transfer_function_value
                , ( bp::arg("i"), bp::arg("segid") )
                , bp::return_self< >() );
        
        }
        { //::SireMol::SegStructureEditor::transferAll
        
            typedef ::SireMol::SegStructureEditor & ( ::SireMol::SegStructureEditor::*transferAll_function_type )( ::SireMol::SegID const & ) ;
            transferAll_function_type transferAll_function_value( &::SireMol::SegStructureEditor::transferAll );
            
            SegStructureEditor_exposer.def( 
                "transferAll"
                , transferAll_function_value
                , ( bp::arg("segid") )
                , bp::return_self< >() );
        
        }
        { //::SireMol::SegStructureEditor::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::SegStructureEditor::typeName );
            
            SegStructureEditor_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMol::SegStructureEditor::what
        
            typedef char const * ( ::SireMol::SegStructureEditor::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireMol::SegStructureEditor::what );
            
            SegStructureEditor_exposer.def( 
                "what"
                , what_function_value );
        
        }
        SegStructureEditor_exposer.staticmethod( "typeName" );
        SegStructureEditor_exposer.def( "__copy__", &__copy__);
        SegStructureEditor_exposer.def( "__deepcopy__", &__copy__);
        SegStructureEditor_exposer.def( "clone", &__copy__);
        SegStructureEditor_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::SegStructureEditor >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SegStructureEditor_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::SegStructureEditor >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SegStructureEditor_exposer.def( "__str__", &__str__< ::SireMol::SegStructureEditor > );
        SegStructureEditor_exposer.def( "__repr__", &__str__< ::SireMol::SegStructureEditor > );
    }

}
