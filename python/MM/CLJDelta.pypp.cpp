// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "CLJDelta.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMol/moleculeview.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "cljdelta.h"

#include "cljdelta.h"

SireMM::CLJDelta __copy__(const SireMM::CLJDelta &other){ return SireMM::CLJDelta(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_CLJDelta_class(){

    { //::SireMM::CLJDelta
        typedef bp::class_< SireMM::CLJDelta > CLJDelta_exposer_t;
        CLJDelta_exposer_t CLJDelta_exposer = CLJDelta_exposer_t( "CLJDelta", bp::init< >() );
        bp::scope CLJDelta_scope( CLJDelta_exposer );
        CLJDelta_exposer.def( bp::init< quint32, SireMM::CLJBoxes const &, QVector< SireMM::CLJBoxIndex > const &, SireMol::MoleculeView const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("idnum"), bp::arg("boxes"), bp::arg("old_atoms"), bp::arg("new_atoms"), bp::arg("map")=SireBase::PropertyMap() )) );
        CLJDelta_exposer.def( bp::init< quint32, SireMM::CLJBoxes const &, QVector< SireMM::CLJBoxIndex > const &, SireMol::MoleculeView const &, SireMM::CLJAtoms::ID_SOURCE, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("idnum"), bp::arg("boxes"), bp::arg("old_atoms"), bp::arg("new_atoms"), bp::arg("source"), bp::arg("map")=SireBase::PropertyMap() )) );
        CLJDelta_exposer.def( bp::init< SireMM::CLJDelta const & >(( bp::arg("other") )) );
        { //::SireMM::CLJDelta::ID
        
            typedef ::quint32 ( ::SireMM::CLJDelta::*ID_function_type )(  ) const;
            ID_function_type ID_function_value( &::SireMM::CLJDelta::ID );
            
            CLJDelta_exposer.def( 
                "ID"
                , ID_function_value );
        
        }
        { //::SireMM::CLJDelta::boxIndex
        
            typedef ::SireMM::CLJBoxIndex const & ( ::SireMM::CLJDelta::*boxIndex_function_type )(  ) const;
            boxIndex_function_type boxIndex_function_value( &::SireMM::CLJDelta::boxIndex );
            
            CLJDelta_exposer.def( 
                "boxIndex"
                , boxIndex_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMM::CLJDelta::boxLength
        
            typedef float ( ::SireMM::CLJDelta::*boxLength_function_type )(  ) const;
            boxLength_function_type boxLength_function_value( &::SireMM::CLJDelta::boxLength );
            
            CLJDelta_exposer.def( 
                "boxLength"
                , boxLength_function_value );
        
        }
        { //::SireMM::CLJDelta::changedAtoms
        
            typedef ::SireMM::CLJAtoms const & ( ::SireMM::CLJDelta::*changedAtoms_function_type )(  ) const;
            changedAtoms_function_type changedAtoms_function_value( &::SireMM::CLJDelta::changedAtoms );
            
            CLJDelta_exposer.def( 
                "changedAtoms"
                , changedAtoms_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMM::CLJDelta::isNull
        
            typedef bool ( ::SireMM::CLJDelta::*isNull_function_type )(  ) const;
            isNull_function_type isNull_function_value( &::SireMM::CLJDelta::isNull );
            
            CLJDelta_exposer.def( 
                "isNull"
                , isNull_function_value );
        
        }
        { //::SireMM::CLJDelta::isSingleBox
        
            typedef bool ( ::SireMM::CLJDelta::*isSingleBox_function_type )(  ) const;
            isSingleBox_function_type isSingleBox_function_value( &::SireMM::CLJDelta::isSingleBox );
            
            CLJDelta_exposer.def( 
                "isSingleBox"
                , isSingleBox_function_value );
        
        }
        { //::SireMM::CLJDelta::maxBox
        
            typedef ::SireMM::CLJBoxIndex ( ::SireMM::CLJDelta::*maxBox_function_type )(  ) const;
            maxBox_function_type maxBox_function_value( &::SireMM::CLJDelta::maxBox );
            
            CLJDelta_exposer.def( 
                "maxBox"
                , maxBox_function_value );
        
        }
        { //::SireMM::CLJDelta::merge
        
            typedef ::SireMM::CLJDelta ( *merge_function_type )( ::SireMM::CLJDelta const *,int,bool );
            merge_function_type merge_function_value( &::SireMM::CLJDelta::merge );
            
            CLJDelta_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("deltas"), bp::arg("count"), bp::arg("changes_only")=(bool)(false) ) );
        
        }
        { //::SireMM::CLJDelta::minBox
        
            typedef ::SireMM::CLJBoxIndex ( ::SireMM::CLJDelta::*minBox_function_type )(  ) const;
            minBox_function_type minBox_function_value( &::SireMM::CLJDelta::minBox );
            
            CLJDelta_exposer.def( 
                "minBox"
                , minBox_function_value );
        
        }
        { //::SireMM::CLJDelta::nBoxX
        
            typedef ::quint8 ( ::SireMM::CLJDelta::*nBoxX_function_type )(  ) const;
            nBoxX_function_type nBoxX_function_value( &::SireMM::CLJDelta::nBoxX );
            
            CLJDelta_exposer.def( 
                "nBoxX"
                , nBoxX_function_value );
        
        }
        { //::SireMM::CLJDelta::nBoxY
        
            typedef ::quint8 ( ::SireMM::CLJDelta::*nBoxY_function_type )(  ) const;
            nBoxY_function_type nBoxY_function_value( &::SireMM::CLJDelta::nBoxY );
            
            CLJDelta_exposer.def( 
                "nBoxY"
                , nBoxY_function_value );
        
        }
        { //::SireMM::CLJDelta::nBoxZ
        
            typedef ::quint8 ( ::SireMM::CLJDelta::*nBoxZ_function_type )(  ) const;
            nBoxZ_function_type nBoxZ_function_value( &::SireMM::CLJDelta::nBoxZ );
            
            CLJDelta_exposer.def( 
                "nBoxZ"
                , nBoxZ_function_value );
        
        }
        { //::SireMM::CLJDelta::nBoxes
        
            typedef int ( ::SireMM::CLJDelta::*nBoxes_function_type )(  ) const;
            nBoxes_function_type nBoxes_function_value( &::SireMM::CLJDelta::nBoxes );
            
            CLJDelta_exposer.def( 
                "nBoxes"
                , nBoxes_function_value );
        
        }
        { //::SireMM::CLJDelta::newAtoms
        
            typedef ::SireMM::CLJAtoms const & ( ::SireMM::CLJDelta::*newAtoms_function_type )(  ) const;
            newAtoms_function_type newAtoms_function_value( &::SireMM::CLJDelta::newAtoms );
            
            CLJDelta_exposer.def( 
                "newAtoms"
                , newAtoms_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMM::CLJDelta::oldAtoms
        
            typedef ::SireMM::CLJAtoms ( ::SireMM::CLJDelta::*oldAtoms_function_type )(  ) const;
            oldAtoms_function_type oldAtoms_function_value( &::SireMM::CLJDelta::oldAtoms );
            
            CLJDelta_exposer.def( 
                "oldAtoms"
                , oldAtoms_function_value );
        
        }
        { //::SireMM::CLJDelta::oldIndicies
        
            typedef ::QVector< SireMM::CLJBoxIndex > const & ( ::SireMM::CLJDelta::*oldIndicies_function_type )(  ) const;
            oldIndicies_function_type oldIndicies_function_value( &::SireMM::CLJDelta::oldIndicies );
            
            CLJDelta_exposer.def( 
                "oldIndicies"
                , oldIndicies_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        CLJDelta_exposer.def( bp::self != bp::self );
        { //::SireMM::CLJDelta::operator=
        
            typedef ::SireMM::CLJDelta & ( ::SireMM::CLJDelta::*assign_function_type )( ::SireMM::CLJDelta const & ) ;
            assign_function_type assign_function_value( &::SireMM::CLJDelta::operator= );
            
            CLJDelta_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        CLJDelta_exposer.def( bp::self == bp::self );
        { //::SireMM::CLJDelta::reconstruct
        
            typedef void ( ::SireMM::CLJDelta::*reconstruct_function_type )( ::quint32,::SireMM::CLJBoxes const &,::QVector< SireMM::CLJBoxIndex > const &,::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) ;
            reconstruct_function_type reconstruct_function_value( &::SireMM::CLJDelta::reconstruct );
            
            CLJDelta_exposer.def( 
                "reconstruct"
                , reconstruct_function_value
                , ( bp::arg("idnum"), bp::arg("boxes"), bp::arg("old_atoms"), bp::arg("new_atoms"), bp::arg("map")=SireBase::PropertyMap() ) );
        
        }
        { //::SireMM::CLJDelta::reconstruct
        
            typedef void ( ::SireMM::CLJDelta::*reconstruct_function_type )( ::quint32,::SireMM::CLJBoxes const &,::QVector< SireMM::CLJBoxIndex > const &,::SireMol::MoleculeView const &,::SireMM::CLJAtoms::ID_SOURCE,::SireBase::PropertyMap const & ) ;
            reconstruct_function_type reconstruct_function_value( &::SireMM::CLJDelta::reconstruct );
            
            CLJDelta_exposer.def( 
                "reconstruct"
                , reconstruct_function_value
                , ( bp::arg("idnum"), bp::arg("boxes"), bp::arg("old_atoms"), bp::arg("new_atoms"), bp::arg("source"), bp::arg("map")=SireBase::PropertyMap() ) );
        
        }
        { //::SireMM::CLJDelta::toString
        
            typedef ::QString ( ::SireMM::CLJDelta::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMM::CLJDelta::toString );
            
            CLJDelta_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMM::CLJDelta::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::CLJDelta::typeName );
            
            CLJDelta_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMM::CLJDelta::what
        
            typedef char const * ( ::SireMM::CLJDelta::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireMM::CLJDelta::what );
            
            CLJDelta_exposer.def( 
                "what"
                , what_function_value );
        
        }
        CLJDelta_exposer.staticmethod( "merge" );
        CLJDelta_exposer.staticmethod( "typeName" );
        CLJDelta_exposer.def( "__copy__", &__copy__);
        CLJDelta_exposer.def( "__deepcopy__", &__copy__);
        CLJDelta_exposer.def( "clone", &__copy__);
        CLJDelta_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::CLJDelta >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CLJDelta_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::CLJDelta >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CLJDelta_exposer.def( "__str__", &__str__< ::SireMM::CLJDelta > );
        CLJDelta_exposer.def( "__repr__", &__str__< ::SireMM::CLJDelta > );
    }

}
