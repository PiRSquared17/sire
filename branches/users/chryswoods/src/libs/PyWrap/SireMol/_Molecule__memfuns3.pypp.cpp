// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "_Molecule__memfuns3.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"

namespace bp = boost::python;

void register_Molecule_memfuns3( Molecule_exposer_t& Molecule_exposer ){

    { //::SireMol::Molecule::coordinates
    
        typedef ::QHash<SireMol::CGAtomID,SireMaths::Vector> ( ::SireMol::Molecule::*coordinates_function_type )( ::QSet<SireMol::CGAtomID> const & ) const;
        
        Molecule_exposer.def( 
            "coordinates"
            , coordinates_function_type( &::SireMol::Molecule::coordinates )
            , ( bp::arg("cgatomids") ) );
    
    }

    { //::SireMol::Molecule::coordinates
    
        typedef ::QHash<SireMol::ResNumAtomID,SireMaths::Vector> ( ::SireMol::Molecule::*coordinates_function_type )( ::QSet<SireMol::ResNumAtomID> const & ) const;
        
        Molecule_exposer.def( 
            "coordinates"
            , coordinates_function_type( &::SireMol::Molecule::coordinates )
            , ( bp::arg("resatomids") ) );
    
    }

    { //::SireMol::Molecule::coordinates
    
        typedef ::QHash<SireMol::ResIDAtomID,SireMaths::Vector> ( ::SireMol::Molecule::*coordinates_function_type )( ::QSet<SireMol::ResIDAtomID> const & ) const;
        
        Molecule_exposer.def( 
            "coordinates"
            , coordinates_function_type( &::SireMol::Molecule::coordinates )
            , ( bp::arg("resatomids") ) );
    
    }

    { //::SireMol::Molecule::coordinates
    
        typedef ::QHash<SireMol::AtomIndex,SireMaths::Vector> ( ::SireMol::Molecule::*coordinates_function_type )( ::QSet<SireMol::AtomIndex> const & ) const;
        
        Molecule_exposer.def( 
            "coordinates"
            , coordinates_function_type( &::SireMol::Molecule::coordinates )
            , ( bp::arg("atoms") ) );
    
    }

    { //::SireMol::Molecule::coordinates
    
        typedef ::QVector<SireMaths::Vector> ( ::SireMol::Molecule::*coordinates_function_type )( ::SireMol::CutGroupID ) ;
        
        Molecule_exposer.def( 
            "coordinates"
            , coordinates_function_type( &::SireMol::Molecule::coordinates )
            , ( bp::arg("cgid") ) );
    
    }

    { //::SireMol::Molecule::coordinates
    
        typedef ::QHash<SireMol::CutGroupID,QVector<SireMaths::Vector> > ( ::SireMol::Molecule::*coordinates_function_type )( ::QSet<SireMol::CutGroupID> const & ) const;
        
        Molecule_exposer.def( 
            "coordinates"
            , coordinates_function_type( &::SireMol::Molecule::coordinates )
            , ( bp::arg("cgids") ) );
    
    }

    { //::SireMol::Molecule::coordinates
    
        typedef ::QVector<SireMaths::Vector> ( ::SireMol::Molecule::*coordinates_function_type )( ::SireMol::ResNum ) ;
        
        Molecule_exposer.def( 
            "coordinates"
            , coordinates_function_type( &::SireMol::Molecule::coordinates )
            , ( bp::arg("resnum") ) );
    
    }

    { //::SireMol::Molecule::coordinates
    
        typedef ::QHash<SireMol::ResNum,QVector<SireMaths::Vector> > ( ::SireMol::Molecule::*coordinates_function_type )( ::QSet<SireMol::ResNum> const & ) const;
        
        Molecule_exposer.def( 
            "coordinates"
            , coordinates_function_type( &::SireMol::Molecule::coordinates )
            , ( bp::arg("resnums") ) );
    
    }

    { //::SireMol::Molecule::coordinates
    
        typedef ::QVector<SireMaths::Vector> ( ::SireMol::Molecule::*coordinates_function_type )( ::SireMol::ResID ) ;
        
        Molecule_exposer.def( 
            "coordinates"
            , coordinates_function_type( &::SireMol::Molecule::coordinates )
            , ( bp::arg("resid") ) );
    
    }

    { //::SireMol::Molecule::coordinates
    
        typedef ::QHash<SireMol::ResID,QVector<SireMaths::Vector> > ( ::SireMol::Molecule::*coordinates_function_type )( ::QSet<SireMol::ResID> const & ) const;
        
        Molecule_exposer.def( 
            "coordinates"
            , coordinates_function_type( &::SireMol::Molecule::coordinates )
            , ( bp::arg("resids") ) );
    
    }

    { //::SireMol::Molecule::cutGroup
    
        typedef ::SireMol::CutGroup ( ::SireMol::Molecule::*cutGroup_function_type )( ::SireMol::CutGroupID ) const;
        
        Molecule_exposer.def( 
            "cutGroup"
            , cutGroup_function_type( &::SireMol::Molecule::cutGroup )
            , ( bp::arg("id") ) );
    
    }

    { //::SireMol::Molecule::cutGroups
    
        typedef ::QVector<SireMol::CutGroup> ( ::SireMol::Molecule::*cutGroups_function_type )(  ) const;
        
        Molecule_exposer.def( 
            "cutGroups"
            , cutGroups_function_type( &::SireMol::Molecule::cutGroups ) );
    
    }

    { //::SireMol::Molecule::cutGroups
    
        typedef ::QHash<SireMol::CutGroupID,SireMol::CutGroup> ( ::SireMol::Molecule::*cutGroups_function_type )( ::SireMol::ResNum ) const;
        
        Molecule_exposer.def( 
            "cutGroups"
            , cutGroups_function_type( &::SireMol::Molecule::cutGroups )
            , ( bp::arg("resnum") ) );
    
    }

    { //::SireMol::Molecule::cutGroups
    
        typedef ::QHash<SireMol::CutGroupID,SireMol::CutGroup> ( ::SireMol::Molecule::*cutGroups_function_type )( ::SireMol::ResID ) const;
        
        Molecule_exposer.def( 
            "cutGroups"
            , cutGroups_function_type( &::SireMol::Molecule::cutGroups )
            , ( bp::arg("resid") ) );
    
    }

    { //::SireMol::Molecule::dihedral
    
        typedef ::SireMaths::Torsion ( ::SireMol::Molecule::*dihedral_function_type )( ::SireMol::Dihedral const & ) const;
        
        Molecule_exposer.def( 
            "dihedral"
            , dihedral_function_type( &::SireMol::Molecule::dihedral )
            , ( bp::arg("dih") ) );
    
    }

    { //::SireMol::Molecule::edit
    
        typedef ::SireMol::EditMol ( ::SireMol::Molecule::*edit_function_type )(  ) const;
        
        Molecule_exposer.def( 
            "edit"
            , edit_function_type( &::SireMol::Molecule::edit ) );
    
    }

    { //::SireMol::Molecule::getProperty
    
        typedef ::SireMol::Property const & ( ::SireMol::Molecule::*getProperty_function_type )( ::QString const & ) const;
        
        Molecule_exposer.def( 
            "getProperty"
            , getProperty_function_type( &::SireMol::Molecule::getProperty )
            , ( bp::arg("name") )
            , bp::return_value_policy< bp::copy_const_reference >() );
    
    }

    { //::SireMol::Molecule::getWeight
    
        typedef double ( ::SireMol::Molecule::*getWeight_function_type )( ::SireMol::AtomIDGroup const &,::SireMol::AtomIDGroup const &,::SireMol::WeightFunction const & ) const;
        
        Molecule_exposer.def( 
            "getWeight"
            , getWeight_function_type( &::SireMol::Molecule::getWeight )
            , ( bp::arg("group0"), bp::arg("group1"), bp::arg("weightfunc") ) );
    
    }

    { //::SireMol::Molecule::idString
    
        typedef ::QString ( ::SireMol::Molecule::*idString_function_type )(  ) const;
        
        Molecule_exposer.def( 
            "idString"
            , idString_function_type( &::SireMol::Molecule::idString ) );
    
    }

    { //::SireMol::Molecule::improper
    
        typedef ::SireMaths::Torsion ( ::SireMol::Molecule::*improper_function_type )( ::SireMol::Improper const & ) const;
        
        Molecule_exposer.def( 
            "improper"
            , improper_function_type( &::SireMol::Molecule::improper )
            , ( bp::arg("improper") ) );
    
    }

}
