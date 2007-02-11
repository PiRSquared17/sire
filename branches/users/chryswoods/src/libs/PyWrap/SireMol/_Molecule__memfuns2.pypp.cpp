// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "_Molecule__memfuns2.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"

namespace bp = boost::python;

void register_Molecule_memfuns2( Molecule_exposer_t& Molecule_exposer ){

    { //::SireMol::Molecule::connectivity
    
        typedef ::SireMol::MoleculeBonds ( ::SireMol::Molecule::*connectivity_function_type )(  ) const;
        
        Molecule_exposer.def( 
            "connectivity"
            , connectivity_function_type( &::SireMol::Molecule::connectivity ) );
    
    }

    { //::SireMol::Molecule::connectivity
    
        typedef ::SireMol::ResidueBonds ( ::SireMol::Molecule::*connectivity_function_type )( ::SireMol::ResNum ) const;
        
        Molecule_exposer.def( 
            "connectivity"
            , connectivity_function_type( &::SireMol::Molecule::connectivity )
            , ( bp::arg("resnum") ) );
    
    }

    { //::SireMol::Molecule::connectivity
    
        typedef ::SireMol::ResidueBonds ( ::SireMol::Molecule::*connectivity_function_type )( ::SireMol::ResID ) const;
        
        Molecule_exposer.def( 
            "connectivity"
            , connectivity_function_type( &::SireMol::Molecule::connectivity )
            , ( bp::arg("resid") ) );
    
    }

    { //::SireMol::Molecule::contains
    
        typedef bool ( ::SireMol::Molecule::*contains_function_type )( ::SireMol::CutGroupID ) const;
        
        Molecule_exposer.def( 
            "contains"
            , contains_function_type( &::SireMol::Molecule::contains )
            , ( bp::arg("cgid") ) );
    
    }

    { //::SireMol::Molecule::contains
    
        typedef bool ( ::SireMol::Molecule::*contains_function_type )( ::SireMol::ResNum ) const;
        
        Molecule_exposer.def( 
            "contains"
            , contains_function_type( &::SireMol::Molecule::contains )
            , ( bp::arg("resnum") ) );
    
    }

    { //::SireMol::Molecule::contains
    
        typedef bool ( ::SireMol::Molecule::*contains_function_type )( ::SireMol::ResID ) const;
        
        Molecule_exposer.def( 
            "contains"
            , contains_function_type( &::SireMol::Molecule::contains )
            , ( bp::arg("resid") ) );
    
    }

    { //::SireMol::Molecule::contains
    
        typedef bool ( ::SireMol::Molecule::*contains_function_type )( ::SireMol::CGAtomID const & ) const;
        
        Molecule_exposer.def( 
            "contains"
            , contains_function_type( &::SireMol::Molecule::contains )
            , ( bp::arg("cgatomid") ) );
    
    }

    { //::SireMol::Molecule::contains
    
        typedef bool ( ::SireMol::Molecule::*contains_function_type )( ::SireMol::IDMolAtom const & ) const;
        
        Molecule_exposer.def( 
            "contains"
            , contains_function_type( &::SireMol::Molecule::contains )
            , ( bp::arg("atomid") ) );
    
    }

    { //::SireMol::Molecule::contains
    
        typedef bool ( ::SireMol::Molecule::*contains_function_type )( ::SireMol::Bond const & ) const;
        
        Molecule_exposer.def( 
            "contains"
            , contains_function_type( &::SireMol::Molecule::contains )
            , ( bp::arg("bond") ) );
    
    }

    { //::SireMol::Molecule::coordGroup
    
        typedef ::SireVol::CoordGroup ( ::SireMol::Molecule::*coordGroup_function_type )( ::SireMol::CutGroupID ) const;
        
        Molecule_exposer.def( 
            "coordGroup"
            , coordGroup_function_type( &::SireMol::Molecule::coordGroup )
            , ( bp::arg("id") ) );
    
    }

    { //::SireMol::Molecule::coordGroups
    
        typedef ::QVector<SireVol::CoordGroup> ( ::SireMol::Molecule::*coordGroups_function_type )(  ) const;
        
        Molecule_exposer.def( 
            "coordGroups"
            , coordGroups_function_type( &::SireMol::Molecule::coordGroups ) );
    
    }

    { //::SireMol::Molecule::coordGroups
    
        typedef ::QHash<SireMol::CutGroupID,SireVol::CoordGroup> ( ::SireMol::Molecule::*coordGroups_function_type )( ::SireMol::ResNum ) const;
        
        Molecule_exposer.def( 
            "coordGroups"
            , coordGroups_function_type( &::SireMol::Molecule::coordGroups )
            , ( bp::arg("resnum") ) );
    
    }

    { //::SireMol::Molecule::coordGroups
    
        typedef ::QHash<SireMol::CutGroupID,SireVol::CoordGroup> ( ::SireMol::Molecule::*coordGroups_function_type )( ::SireMol::ResID ) const;
        
        Molecule_exposer.def( 
            "coordGroups"
            , coordGroups_function_type( &::SireMol::Molecule::coordGroups )
            , ( bp::arg("resid") ) );
    
    }

    { //::SireMol::Molecule::coordGroups
    
        typedef ::QHash<SireMol::CutGroupID,SireVol::CoordGroup> ( ::SireMol::Molecule::*coordGroups_function_type )( ::QSet<SireMol::CutGroupID> const & ) const;
        
        Molecule_exposer.def( 
            "coordGroups"
            , coordGroups_function_type( &::SireMol::Molecule::coordGroups )
            , ( bp::arg("cgids") ) );
    
    }

    { //::SireMol::Molecule::coordGroups
    
        typedef ::QHash<SireMol::CutGroupID,SireVol::CoordGroup> ( ::SireMol::Molecule::*coordGroups_function_type )( ::QSet<SireMol::ResNum> const & ) const;
        
        Molecule_exposer.def( 
            "coordGroups"
            , coordGroups_function_type( &::SireMol::Molecule::coordGroups )
            , ( bp::arg("resnums") ) );
    
    }

    { //::SireMol::Molecule::coordGroups
    
        typedef ::QHash<SireMol::CutGroupID,SireVol::CoordGroup> ( ::SireMol::Molecule::*coordGroups_function_type )( ::QSet<SireMol::ResID> const & ) const;
        
        Molecule_exposer.def( 
            "coordGroups"
            , coordGroups_function_type( &::SireMol::Molecule::coordGroups )
            , ( bp::arg("resids") ) );
    
    }

    { //::SireMol::Molecule::coordinates
    
        typedef ::SireMaths::Vector ( ::SireMol::Molecule::*coordinates_function_type )( ::SireMol::CGAtomID const & ) const;
        
        Molecule_exposer.def( 
            "coordinates"
            , coordinates_function_type( &::SireMol::Molecule::coordinates )
            , ( bp::arg("cgatomid") ) );
    
    }

    { //::SireMol::Molecule::coordinates
    
        typedef ::SireMaths::Vector ( ::SireMol::Molecule::*coordinates_function_type )( ::SireMol::IDMolAtom const & ) const;
        
        Molecule_exposer.def( 
            "coordinates"
            , coordinates_function_type( &::SireMol::Molecule::coordinates )
            , ( bp::arg("atomid") ) );
    
    }

    { //::SireMol::Molecule::coordinates
    
        typedef ::QVector<SireMaths::Vector> ( ::SireMol::Molecule::*coordinates_function_type )(  ) const;
        
        Molecule_exposer.def( 
            "coordinates"
            , coordinates_function_type( &::SireMol::Molecule::coordinates ) );
    
    }

    { //::SireMol::Molecule::coordinates
    
        typedef ::QHash<SireMol::AtomID,SireMaths::Vector> ( ::SireMol::Molecule::*coordinates_function_type )( ::QSet<SireMol::AtomID> const & ) const;
        
        Molecule_exposer.def( 
            "coordinates"
            , coordinates_function_type( &::SireMol::Molecule::coordinates )
            , ( bp::arg("atomids") ) );
    
    }

}
