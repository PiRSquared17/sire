// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "Residue.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMol/moleculedata.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"
#include "SireVol/space.h"

namespace bp = boost::python;

SireMol::Residue __copy__(const SireMol::Residue &other){ return SireMol::Residue(other); }

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireMol::Residue&){ return "SireMol::Residue";}

void register_Residue_class(){

    bp::class_< SireMol::Residue, bp::bases< SireMol::MoleculeView > >( "Residue" )    
        .def( bp::init< >() )    
        .def( bp::init< SireMol::Molecule const &, SireMol::ResNum >(( bp::arg("mol"), bp::arg("resnum") )) )    
        .def( bp::init< SireMol::Molecule const &, SireMol::ResID >(( bp::arg("mol"), bp::arg("resid") )) )    
        .def( bp::init< SireMol::Molecule const &, QString const & >(( bp::arg("mol"), bp::arg("resname") )) )    
        .def( bp::init< SireMol::NewAtom const & >(( bp::arg("atom") )) )    
        .def( 
            "ID"
            , &::SireMol::Residue::ID )    
        .def( 
            "angle"
            , &::SireMol::Residue::angle
            , ( bp::arg("ang") ) )    
        .def( 
            "at"
            , (::SireMol::Atom ( ::SireMol::Residue::* )( ::SireMol::AtomID ) const)( &::SireMol::Residue::at )
            , ( bp::arg("i") ) )    
        .def( 
            "at"
            , (::SireMol::Atom ( ::SireMol::Residue::* )( ::QString const & ) const)( &::SireMol::Residue::at )
            , ( bp::arg("atomname") ) )    
        .def( 
            "at"
            , (::SireMol::Atom ( ::SireMol::Residue::* )( ::SireMol::AtomIndex const & ) const)( &::SireMol::Residue::at )
            , ( bp::arg("atom") ) )    
        .def( 
            "at"
            , (::SireMol::CutGroup ( ::SireMol::Residue::* )( ::SireMol::CutGroupID ) const)( &::SireMol::Residue::at )
            , ( bp::arg("cgid") ) )    
        .def( 
            "at"
            , (::SireMol::CutGroup ( ::SireMol::Residue::* )( ::SireMol::CutGroupNum ) const)( &::SireMol::Residue::at )
            , ( bp::arg("cgnum") ) )    
        .def( 
            "atom"
            , (::SireMol::Atom ( ::SireMol::Residue::* )( ::SireMol::AtomID ) const)( &::SireMol::Residue::atom )
            , ( bp::arg("i") ) )    
        .def( 
            "atom"
            , (::SireMol::Atom ( ::SireMol::Residue::* )( ::QString const & ) const)( &::SireMol::Residue::atom )
            , ( bp::arg("atomname") ) )    
        .def( 
            "atom"
            , (::SireMol::Atom ( ::SireMol::Residue::* )( ::SireMol::AtomIndex const & ) const)( &::SireMol::Residue::atom )
            , ( bp::arg("atom") ) )    
        .def( 
            "atomNames"
            , &::SireMol::Residue::atomNames )    
        .def( 
            "atoms"
            , (::QVector<SireMol::Atom> ( ::SireMol::Residue::* )(  ) const)( &::SireMol::Residue::atoms ) )    
        .def( 
            "atoms"
            , (::QHash<SireMol::AtomID,SireMol::Atom> ( ::SireMol::Residue::* )( ::QSet<SireMol::AtomID> const & ) const)( &::SireMol::Residue::atoms )
            , ( bp::arg("atomids") ) )    
        .def( 
            "atoms"
            , (::QHash<QString,SireMol::Atom> ( ::SireMol::Residue::* )( ::QSet<QString> const & ) const)( &::SireMol::Residue::atoms )
            , ( bp::arg("atomnames") ) )    
        .def( 
            "atoms"
            , (::QHash<SireMol::AtomIndex,SireMol::Atom> ( ::SireMol::Residue::* )( ::QSet<SireMol::AtomIndex> const & ) const)( &::SireMol::Residue::atoms )
            , ( bp::arg("atms") ) )    
        .def( 
            "bond"
            , &::SireMol::Residue::bond
            , ( bp::arg("bnd") ) )    
        .def( 
            "bondedResidues"
            , &::SireMol::Residue::bondedResidues )    
        .def( 
            "connectivity"
            , &::SireMol::Residue::connectivity )    
        .def( 
            "contains"
            , (bool ( ::SireMol::Residue::* )( ::SireMol::CutGroupID ) const)( &::SireMol::Residue::contains )
            , ( bp::arg("cgid") ) )    
        .def( 
            "contains"
            , (bool ( ::SireMol::Residue::* )( ::SireMol::AtomID ) const)( &::SireMol::Residue::contains )
            , ( bp::arg("atomid") ) )    
        .def( 
            "contains"
            , (bool ( ::SireMol::Residue::* )( ::QString const & ) const)( &::SireMol::Residue::contains )
            , ( bp::arg("atomname") ) )    
        .def( 
            "contains"
            , (bool ( ::SireMol::Residue::* )( ::SireMol::AtomIndex const & ) const)( &::SireMol::Residue::contains )
            , ( bp::arg("atom") ) )    
        .def( 
            "contains"
            , (bool ( ::SireMol::Residue::* )( ::SireMol::Bond const & ) const)( &::SireMol::Residue::contains )
            , ( bp::arg("bond") ) )    
        .def( 
            "coordGroup"
            , &::SireMol::Residue::coordGroup
            , ( bp::arg("cgid") ) )    
        .def( 
            "coordGroups"
            , &::SireMol::Residue::coordGroups )    
        .def( 
            "coordinates"
            , (::SireMaths::Vector ( ::SireMol::Residue::* )( ::SireMol::AtomID ) const)( &::SireMol::Residue::coordinates )
            , ( bp::arg("i") ) )    
        .def( 
            "coordinates"
            , (::SireMaths::Vector ( ::SireMol::Residue::* )( ::QString const & ) const)( &::SireMol::Residue::coordinates )
            , ( bp::arg("atomname") ) )    
        .def( 
            "coordinates"
            , (::SireMaths::Vector ( ::SireMol::Residue::* )( ::SireMol::AtomIndex const & ) const)( &::SireMol::Residue::coordinates )
            , ( bp::arg("atom") ) )    
        .def( 
            "coordinates"
            , (::QVector<SireMaths::Vector> ( ::SireMol::Residue::* )(  ) const)( &::SireMol::Residue::coordinates ) )    
        .def( 
            "coordinates"
            , (::QHash<SireMol::AtomID,SireMaths::Vector> ( ::SireMol::Residue::* )( ::QSet<SireMol::AtomID> const & ) const)( &::SireMol::Residue::coordinates )
            , ( bp::arg("atomids") ) )    
        .def( 
            "coordinates"
            , (::QHash<QString,SireMaths::Vector> ( ::SireMol::Residue::* )( ::QSet<QString> const & ) const)( &::SireMol::Residue::coordinates )
            , ( bp::arg("atomnames") ) )    
        .def( 
            "coordinates"
            , (::QHash<SireMol::AtomIndex,SireMaths::Vector> ( ::SireMol::Residue::* )( ::QSet<SireMol::AtomIndex> const & ) const)( &::SireMol::Residue::coordinates )
            , ( bp::arg("atms") ) )    
        .def( 
            "cutGroup"
            , &::SireMol::Residue::cutGroup
            , ( bp::arg("cgid") ) )    
        .def( 
            "cutGroups"
            , &::SireMol::Residue::cutGroups )    
        .def( 
            "dihedral"
            , &::SireMol::Residue::dihedral
            , ( bp::arg("dih") ) )    
        .def( 
            "getWeight"
            , (double ( ::SireMol::Residue::* )( ::QStringList const &,::QStringList const &,::SireMol::WeightFunction const & ) const)( &::SireMol::Residue::getWeight )
            , ( bp::arg("group0"), bp::arg("group1"), bp::arg("weightfunc") ) )    
        .def( 
            "getWeight"
            , (double ( ::SireMol::Residue::* )( ::QSet<SireMol::AtomIndex> const &,::QSet<SireMol::AtomIndex> const &,::SireMol::WeightFunction const & ) const)( &::SireMol::Residue::getWeight )
            , ( bp::arg("group0"), bp::arg("group1"), bp::arg("weightfunc") ) )    
        .def( 
            "idString"
            , &::SireMol::Residue::idString )    
        .def( 
            "improper"
            , &::SireMol::Residue::improper
            , ( bp::arg("improper") ) )    
        .def( 
            "info"
            , &::SireMol::Residue::info
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "isEmpty"
            , &::SireMol::Residue::isEmpty )    
        .def( 
            "measure"
            , (double ( ::SireMol::Residue::* )( ::SireMol::Bond const & ) const)( &::SireMol::Residue::measure )
            , ( bp::arg("bnd") ) )    
        .def( 
            "measure"
            , (::SireUnits::Dimension::Angle ( ::SireMol::Residue::* )( ::SireMol::Angle const & ) const)( &::SireMol::Residue::measure )
            , ( bp::arg("ang") ) )    
        .def( 
            "measure"
            , (::SireUnits::Dimension::Angle ( ::SireMol::Residue::* )( ::SireMol::Dihedral const & ) const)( &::SireMol::Residue::measure )
            , ( bp::arg("dih") ) )    
        .def( 
            "measure"
            , (::SireUnits::Dimension::Angle ( ::SireMol::Residue::* )( ::SireMol::Improper const & ) const)( &::SireMol::Residue::measure )
            , ( bp::arg("improper") ) )    
        .def( 
            "nAtoms"
            , (int ( ::SireMol::Residue::* )(  ) const)( &::SireMol::Residue::nAtoms ) )    
        .def( 
            "nAtoms"
            , (int ( ::SireMol::Residue::* )( ::SireMol::CutGroupID ) const)( &::SireMol::Residue::nAtoms )
            , ( bp::arg("cgid") ) )    
        .def( 
            "nBonds"
            , &::SireMol::Residue::nBonds )    
        .def( 
            "nCutGroups"
            , &::SireMol::Residue::nCutGroups )    
        .def( 
            "nInterBonds"
            , &::SireMol::Residue::nInterBonds )    
        .def( 
            "nIntraBonds"
            , &::SireMol::Residue::nIntraBonds )    
        .def( 
            "name"
            , &::SireMol::Residue::name )    
        .def( 
            "number"
            , &::SireMol::Residue::number )    
        .def( bp::self != bp::self )    
        .def( bp::self == bp::self )    
        .def( 
            "__getitem__"
            , (::SireMol::Atom ( ::SireMol::Residue::* )( ::SireMol::AtomID ) const)( &::SireMol::Residue::operator[] )
            , ( bp::arg("i") ) )    
        .def( 
            "__getitem__"
            , (::SireMol::Atom ( ::SireMol::Residue::* )( ::QString const & ) const)( &::SireMol::Residue::operator[] )
            , ( bp::arg("atomname") ) )    
        .def( 
            "__getitem__"
            , (::SireMol::Atom ( ::SireMol::Residue::* )( ::SireMol::AtomIndex const & ) const)( &::SireMol::Residue::operator[] )
            , ( bp::arg("atom") ) )    
        .def( 
            "__getitem__"
            , (::SireMol::CutGroup ( ::SireMol::Residue::* )( ::SireMol::CutGroupID ) const)( &::SireMol::Residue::operator[] )
            , ( bp::arg("cgid") ) )    
        .def( 
            "__getitem__"
            , (::SireMol::CutGroup ( ::SireMol::Residue::* )( ::SireMol::CutGroupNum ) const)( &::SireMol::Residue::operator[] )
            , ( bp::arg("cgnum") ) )    
        .def( 
            "resName"
            , &::SireMol::Residue::resName )    
        .def( 
            "resNum"
            , &::SireMol::Residue::resNum )    
        .def( 
            "residuesBondedTo"
            , (::QHash<SireMol::ResNum,SireMol::Residue> ( ::SireMol::Residue::* )( ::SireMol::AtomID ) const)( &::SireMol::Residue::residuesBondedTo )
            , ( bp::arg("atom") ) )    
        .def( 
            "residuesBondedTo"
            , (::QHash<SireMol::ResNum,SireMol::Residue> ( ::SireMol::Residue::* )( ::QString const & ) const)( &::SireMol::Residue::residuesBondedTo )
            , ( bp::arg("atomname") ) )    
        .def( 
            "residuesBondedTo"
            , (::QHash<SireMol::ResNum,SireMol::Residue> ( ::SireMol::Residue::* )( ::SireMol::AtomIndex const & ) const)( &::SireMol::Residue::residuesBondedTo )
            , ( bp::arg("atom") ) )    
        .def( 
            "version"
            , &::SireMol::Residue::version
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMol::Residue >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMol::Residue >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
