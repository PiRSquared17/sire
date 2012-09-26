// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "InternalParameters3D.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireFF/errors.h"

#include "SireMol/cgidx.h"

#include "SireMol/molecule.h"

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "SireVol/coordgroup.h"

#include "internalparameters.h"

#include "sireglobal.h"

#include "tostring.h"

#include <algorithm>

#include "internalparameters.h"

SireMM::InternalParameters3D __copy__(const SireMM::InternalParameters3D &other){ return SireMM::InternalParameters3D(other); }

#include "Qt/qdatastream.hpp"

const char* pvt_get_name(const SireMM::InternalParameters3D&){ return "SireMM::InternalParameters3D";}

void register_InternalParameters3D_class(){

    { //::SireMM::InternalParameters3D
        typedef bp::class_< SireMM::InternalParameters3D, bp::bases< SireMM::InternalParameters > > InternalParameters3D_exposer_t;
        InternalParameters3D_exposer_t InternalParameters3D_exposer = InternalParameters3D_exposer_t( "InternalParameters3D", bp::init< >() );
        bp::scope InternalParameters3D_scope( InternalParameters3D_exposer );
        InternalParameters3D_exposer.def( bp::init< SireMol::PartialMolecule const &, SireBase::PropertyName const &, SireBase::PropertyName const &, SireBase::PropertyName const &, SireBase::PropertyName const &, SireBase::PropertyName const &, SireBase::PropertyName const &, SireBase::PropertyName const &, SireBase::PropertyName const &, SireBase::PropertyName const &, SireBase::PropertyName const &, bool >(( bp::arg("molecule"), bp::arg("coords_property"), bp::arg("bond_params"), bp::arg("angle_params"), bp::arg("dihedral_params"), bp::arg("improper_params"), bp::arg("ub_params"), bp::arg("ss_params"), bp::arg("sb_params"), bp::arg("bb_params"), bp::arg("sbt_params"), bp::arg("isstrict") )) );
        InternalParameters3D_exposer.def( bp::init< SireMM::InternalParameters3D const & >(( bp::arg("other") )) );
        { //::SireMM::InternalParameters3D::addChangedGroups
        
            typedef void ( ::SireMM::InternalParameters3D::*addChangedGroups_function_type )( ::SireMM::InternalParameters3D const &,::QSet< unsigned int > & ) const;
            addChangedGroups_function_type addChangedGroups_function_value( &::SireMM::InternalParameters3D::addChangedGroups );
            
            InternalParameters3D_exposer.def( 
                "addChangedGroups"
                , addChangedGroups_function_value
                , ( bp::arg("other"), bp::arg("changed_groups") ) );
        
        }
        { //::SireMM::InternalParameters3D::applyMask
        
            typedef ::SireMM::InternalParameters3D ( ::SireMM::InternalParameters3D::*applyMask_function_type )( ::QSet< unsigned int > const & ) const;
            applyMask_function_type applyMask_function_value( &::SireMM::InternalParameters3D::applyMask );
            
            InternalParameters3D_exposer.def( 
                "applyMask"
                , applyMask_function_value
                , ( bp::arg("cgidxs") ) );
        
        }
        { //::SireMM::InternalParameters3D::atomicCoordinates
        
            typedef ::SireVol::CoordGroupArray const & ( ::SireMM::InternalParameters3D::*atomicCoordinates_function_type )(  ) const;
            atomicCoordinates_function_type atomicCoordinates_function_value( &::SireMM::InternalParameters3D::atomicCoordinates );
            
            InternalParameters3D_exposer.def( 
                "atomicCoordinates"
                , atomicCoordinates_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMM::InternalParameters3D::changedAllGroups
        
            typedef bool ( ::SireMM::InternalParameters3D::*changedAllGroups_function_type )( ::SireMM::InternalParameters3D const & ) const;
            changedAllGroups_function_type changedAllGroups_function_value( &::SireMM::InternalParameters3D::changedAllGroups );
            
            InternalParameters3D_exposer.def( 
                "changedAllGroups"
                , changedAllGroups_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireMM::InternalParameters3D::getChangedGroups
        
            typedef ::QSet< unsigned int > ( ::SireMM::InternalParameters3D::*getChangedGroups_function_type )( ::SireMM::InternalParameters3D const & ) const;
            getChangedGroups_function_type getChangedGroups_function_value( &::SireMM::InternalParameters3D::getChangedGroups );
            
            InternalParameters3D_exposer.def( 
                "getChangedGroups"
                , getChangedGroups_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireMM::InternalParameters3D::nCutGroups
        
            typedef int ( ::SireMM::InternalParameters3D::*nCutGroups_function_type )(  ) const;
            nCutGroups_function_type nCutGroups_function_value( &::SireMM::InternalParameters3D::nCutGroups );
            
            InternalParameters3D_exposer.def( 
                "nCutGroups"
                , nCutGroups_function_value );
        
        }
        InternalParameters3D_exposer.def( bp::self != bp::self );
        { //::SireMM::InternalParameters3D::operator=
        
            typedef ::SireMM::InternalParameters3D & ( ::SireMM::InternalParameters3D::*assign_function_type )( ::SireMM::InternalParameters3D const & ) ;
            assign_function_type assign_function_value( &::SireMM::InternalParameters3D::operator= );
            
            InternalParameters3D_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        InternalParameters3D_exposer.def( bp::self == bp::self );
        { //::SireMM::InternalParameters3D::setAtomicCoordinates
        
            typedef void ( ::SireMM::InternalParameters3D::*setAtomicCoordinates_function_type )( ::SireFF::detail::AtomicCoords3D const & ) ;
            setAtomicCoordinates_function_type setAtomicCoordinates_function_value( &::SireMM::InternalParameters3D::setAtomicCoordinates );
            
            InternalParameters3D_exposer.def( 
                "setAtomicCoordinates"
                , setAtomicCoordinates_function_value
                , ( bp::arg("coords") ) );
        
        }
        { //::SireMM::InternalParameters3D::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::InternalParameters3D::typeName );
            
            InternalParameters3D_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMM::InternalParameters3D::what
        
            typedef char const * ( ::SireMM::InternalParameters3D::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireMM::InternalParameters3D::what );
            
            InternalParameters3D_exposer.def( 
                "what"
                , what_function_value );
        
        }
        InternalParameters3D_exposer.staticmethod( "typeName" );
        InternalParameters3D_exposer.def( "__copy__", &__copy__);
        InternalParameters3D_exposer.def( "__deepcopy__", &__copy__);
        InternalParameters3D_exposer.def( "clone", &__copy__);
        InternalParameters3D_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::InternalParameters3D >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        InternalParameters3D_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::InternalParameters3D >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        InternalParameters3D_exposer.def( "__str__", &pvt_get_name);
        InternalParameters3D_exposer.def( "__repr__", &pvt_get_name);
    }

}
