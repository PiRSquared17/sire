// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "AtomicVelocityWorkspace.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/quickcopy.hpp"

#include "SireMol/atomcoords.h"

#include "SireMol/errors.h"

#include "SireMol/molecule.h"

#include "SireMol/moleculeview.h"

#include "SireMol/moleditor.h"

#include "SireMol/partialmolecule.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "integratorworkspace.h"

#include "velocitygenerator.h"

#include "integratorworkspace.h"

SireMove::AtomicVelocityWorkspace __copy__(const SireMove::AtomicVelocityWorkspace &other){ return SireMove::AtomicVelocityWorkspace(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_AtomicVelocityWorkspace_class(){

    { //::SireMove::AtomicVelocityWorkspace
        typedef bp::class_< SireMove::AtomicVelocityWorkspace, bp::bases< SireMove::IntegratorWorkspace, SireBase::Property > > AtomicVelocityWorkspace_exposer_t;
        AtomicVelocityWorkspace_exposer_t AtomicVelocityWorkspace_exposer = AtomicVelocityWorkspace_exposer_t( "AtomicVelocityWorkspace" );
        bp::scope AtomicVelocityWorkspace_scope( AtomicVelocityWorkspace_exposer );
        AtomicVelocityWorkspace_exposer.def( bp::init< >() );
        AtomicVelocityWorkspace_exposer.def( bp::init< SireMol::MoleculeGroup const & >(( bp::arg("molgroup") )) );
        AtomicVelocityWorkspace_exposer.def( bp::init< SireMove::AtomicVelocityWorkspace const & >(( bp::arg("other") )) );
        { //::SireMove::AtomicVelocityWorkspace::coordinateArray
        
            typedef ::QVector< SireMaths::Vector > const & ( ::SireMove::AtomicVelocityWorkspace::*coordinateArray_function_type )( int ) const;
            coordinateArray_function_type coordinateArray_function_value( &::SireMove::AtomicVelocityWorkspace::coordinateArray );
            
            AtomicVelocityWorkspace_exposer.def( 
                "coordinateArray"
                , coordinateArray_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMove::AtomicVelocityWorkspace::count
        
            typedef int ( ::SireMove::AtomicVelocityWorkspace::*count_function_type )(  ) const;
            count_function_type count_function_value( &::SireMove::AtomicVelocityWorkspace::count );
            
            AtomicVelocityWorkspace_exposer.def( 
                "count"
                , count_function_value );
        
        }
        { //::SireMove::AtomicVelocityWorkspace::forceArray
        
            typedef ::QVector< SireMaths::Vector > const & ( ::SireMove::AtomicVelocityWorkspace::*forceArray_function_type )( int ) const;
            forceArray_function_type forceArray_function_value( &::SireMove::AtomicVelocityWorkspace::forceArray );
            
            AtomicVelocityWorkspace_exposer.def( 
                "forceArray"
                , forceArray_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMove::AtomicVelocityWorkspace::forces
        
            typedef ::QHash< SireMol::MolNum, SireMol::AtomProperty< SireBase::Vector3D< SireUnits::Dimension::PhysUnit< 1, 1, -0x00000000000000002, 0, 0, 0, 0 > > > > ( ::SireMove::AtomicVelocityWorkspace::*forces_function_type )(  ) const;
            forces_function_type forces_function_value( &::SireMove::AtomicVelocityWorkspace::forces );
            
            AtomicVelocityWorkspace_exposer.def( 
                "forces"
                , forces_function_value );
        
        }
        { //::SireMove::AtomicVelocityWorkspace::forces
        
            typedef ::SireMol::AtomForces ( ::SireMove::AtomicVelocityWorkspace::*forces_function_type )( ::SireMol::MolID const & ) const;
            forces_function_type forces_function_value( &::SireMove::AtomicVelocityWorkspace::forces );
            
            AtomicVelocityWorkspace_exposer.def( 
                "forces"
                , forces_function_value
                , ( bp::arg("molid") ) );
        
        }
        { //::SireMove::AtomicVelocityWorkspace::kineticEnergy
        
            typedef ::SireUnits::Dimension::MolarEnergy ( ::SireMove::AtomicVelocityWorkspace::*kineticEnergy_function_type )(  ) const;
            kineticEnergy_function_type kineticEnergy_function_value( &::SireMove::AtomicVelocityWorkspace::kineticEnergy );
            
            AtomicVelocityWorkspace_exposer.def( 
                "kineticEnergy"
                , kineticEnergy_function_value );
        
        }
        { //::SireMove::AtomicVelocityWorkspace::kineticEnergy
        
            typedef ::SireUnits::Dimension::MolarEnergy ( ::SireMove::AtomicVelocityWorkspace::*kineticEnergy_function_type )( ::SireMol::MoleculeView const & ) const;
            kineticEnergy_function_type kineticEnergy_function_value( &::SireMove::AtomicVelocityWorkspace::kineticEnergy );
            
            AtomicVelocityWorkspace_exposer.def( 
                "kineticEnergy"
                , kineticEnergy_function_value
                , ( bp::arg("molview") ) );
        
        }
        { //::SireMove::AtomicVelocityWorkspace::masses
        
            typedef ::QHash< SireMol::MolNum, SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 0, 0, 0, 0, -0x00000000000000001, 0 > > > ( ::SireMove::AtomicVelocityWorkspace::*masses_function_type )(  ) const;
            masses_function_type masses_function_value( &::SireMove::AtomicVelocityWorkspace::masses );
            
            AtomicVelocityWorkspace_exposer.def( 
                "masses"
                , masses_function_value );
        
        }
        { //::SireMove::AtomicVelocityWorkspace::masses
        
            typedef ::SireMol::AtomMasses ( ::SireMove::AtomicVelocityWorkspace::*masses_function_type )( ::SireMol::MolID const & ) const;
            masses_function_type masses_function_value( &::SireMove::AtomicVelocityWorkspace::masses );
            
            AtomicVelocityWorkspace_exposer.def( 
                "masses"
                , masses_function_value
                , ( bp::arg("molid") ) );
        
        }
        { //::SireMove::AtomicVelocityWorkspace::nMolecules
        
            typedef int ( ::SireMove::AtomicVelocityWorkspace::*nMolecules_function_type )(  ) const;
            nMolecules_function_type nMolecules_function_value( &::SireMove::AtomicVelocityWorkspace::nMolecules );
            
            AtomicVelocityWorkspace_exposer.def( 
                "nMolecules"
                , nMolecules_function_value );
        
        }
        AtomicVelocityWorkspace_exposer.def( bp::self != bp::self );
        { //::SireMove::AtomicVelocityWorkspace::operator=
        
            typedef ::SireMove::AtomicVelocityWorkspace & ( ::SireMove::AtomicVelocityWorkspace::*assign_function_type )( ::SireMove::AtomicVelocityWorkspace const & ) ;
            assign_function_type assign_function_value( &::SireMove::AtomicVelocityWorkspace::operator= );
            
            AtomicVelocityWorkspace_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        AtomicVelocityWorkspace_exposer.def( bp::self == bp::self );
        { //::SireMove::AtomicVelocityWorkspace::reciprocalMassArray
        
            typedef ::QVector< double > const & ( ::SireMove::AtomicVelocityWorkspace::*reciprocalMassArray_function_type )( int ) const;
            reciprocalMassArray_function_type reciprocalMassArray_function_value( &::SireMove::AtomicVelocityWorkspace::reciprocalMassArray );
            
            AtomicVelocityWorkspace_exposer.def( 
                "reciprocalMassArray"
                , reciprocalMassArray_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMove::AtomicVelocityWorkspace::setCoordinates
        
            typedef void ( ::SireMove::AtomicVelocityWorkspace::*setCoordinates_function_type )( int,::QVector< SireMaths::Vector > const & ) ;
            setCoordinates_function_type setCoordinates_function_value( &::SireMove::AtomicVelocityWorkspace::setCoordinates );
            
            AtomicVelocityWorkspace_exposer.def( 
                "setCoordinates"
                , setCoordinates_function_value
                , ( bp::arg("i"), bp::arg("new_coords") ) );
        
        }
        { //::SireMove::AtomicVelocityWorkspace::setVelocities
        
            typedef void ( ::SireMove::AtomicVelocityWorkspace::*setVelocities_function_type )( int,::QVector< SireMaths::Vector > const & ) ;
            setVelocities_function_type setVelocities_function_value( &::SireMove::AtomicVelocityWorkspace::setVelocities );
            
            AtomicVelocityWorkspace_exposer.def( 
                "setVelocities"
                , setVelocities_function_value
                , ( bp::arg("i"), bp::arg("new_velocities") ) );
        
        }
        { //::SireMove::AtomicVelocityWorkspace::size
        
            typedef int ( ::SireMove::AtomicVelocityWorkspace::*size_function_type )(  ) const;
            size_function_type size_function_value( &::SireMove::AtomicVelocityWorkspace::size );
            
            AtomicVelocityWorkspace_exposer.def( 
                "size"
                , size_function_value );
        
        }
        { //::SireMove::AtomicVelocityWorkspace::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::AtomicVelocityWorkspace::typeName );
            
            AtomicVelocityWorkspace_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMove::AtomicVelocityWorkspace::updateFrom
        
            typedef void ( ::SireMove::AtomicVelocityWorkspace::*updateFrom_function_type )( ::SireSystem::System &,::SireCAS::Symbol const &,::SireMove::VelocityGenerator const &,::SireBase::PropertyMap const & ) ;
            updateFrom_function_type updateFrom_function_value( &::SireMove::AtomicVelocityWorkspace::updateFrom );
            
            AtomicVelocityWorkspace_exposer.def( 
                "updateFrom"
                , updateFrom_function_value
                , ( bp::arg("system"), bp::arg("nrg_component"), bp::arg("velgen"), bp::arg("map")=SireBase::PropertyMap() ) );
        
        }
        { //::SireMove::AtomicVelocityWorkspace::updateSystem
        
            typedef void ( ::SireMove::AtomicVelocityWorkspace::*updateSystem_function_type )( ::SireSystem::System &,::SireCAS::Symbol const &,::SireBase::PropertyMap const & ) ;
            updateSystem_function_type updateSystem_function_value( &::SireMove::AtomicVelocityWorkspace::updateSystem );
            
            AtomicVelocityWorkspace_exposer.def( 
                "updateSystem"
                , updateSystem_function_value
                , ( bp::arg("system"), bp::arg("nrg_component"), bp::arg("map")=SireBase::PropertyMap() ) );
        
        }
        { //::SireMove::AtomicVelocityWorkspace::velocities
        
            typedef ::QHash< SireMol::MolNum, SireMol::AtomProperty< SireBase::Vector3D< SireUnits::Dimension::PhysUnit< 0, 1, -0x00000000000000001, 0, 0, 0, 0 > > > > ( ::SireMove::AtomicVelocityWorkspace::*velocities_function_type )(  ) const;
            velocities_function_type velocities_function_value( &::SireMove::AtomicVelocityWorkspace::velocities );
            
            AtomicVelocityWorkspace_exposer.def( 
                "velocities"
                , velocities_function_value );
        
        }
        { //::SireMove::AtomicVelocityWorkspace::velocities
        
            typedef ::SireMol::AtomVelocities ( ::SireMove::AtomicVelocityWorkspace::*velocities_function_type )( ::SireMol::MolID const & ) const;
            velocities_function_type velocities_function_value( &::SireMove::AtomicVelocityWorkspace::velocities );
            
            AtomicVelocityWorkspace_exposer.def( 
                "velocities"
                , velocities_function_value
                , ( bp::arg("molid") ) );
        
        }
        { //::SireMove::AtomicVelocityWorkspace::velocityArray
        
            typedef ::QVector< SireMaths::Vector > const & ( ::SireMove::AtomicVelocityWorkspace::*velocityArray_function_type )( int ) const;
            velocityArray_function_type velocityArray_function_value( &::SireMove::AtomicVelocityWorkspace::velocityArray );
            
            AtomicVelocityWorkspace_exposer.def( 
                "velocityArray"
                , velocityArray_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        AtomicVelocityWorkspace_exposer.staticmethod( "typeName" );
        AtomicVelocityWorkspace_exposer.def( "__copy__", &__copy__);
        AtomicVelocityWorkspace_exposer.def( "__deepcopy__", &__copy__);
        AtomicVelocityWorkspace_exposer.def( "clone", &__copy__);
        AtomicVelocityWorkspace_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::AtomicVelocityWorkspace >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomicVelocityWorkspace_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::AtomicVelocityWorkspace >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomicVelocityWorkspace_exposer.def( "__str__", &__str__< ::SireMove::AtomicVelocityWorkspace > );
        AtomicVelocityWorkspace_exposer.def( "__repr__", &__str__< ::SireMove::AtomicVelocityWorkspace > );
    }

}
