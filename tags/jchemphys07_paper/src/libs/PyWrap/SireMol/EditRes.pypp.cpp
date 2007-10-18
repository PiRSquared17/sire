// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "EditRes.pypp.hpp"
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

SireMol::EditRes __copy__(const SireMol::EditRes &other){ return SireMol::EditRes(other); }

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireMol::EditRes&){ return "SireMol::EditRes";}

void register_EditRes_class(){

    bp::class_< SireMol::EditRes >( "EditRes" )    
        .def( bp::init< >() )    
        .def( bp::init< SireMol::EditMol const &, SireMol::ResNum >(( bp::arg("editmol"), bp::arg("resnum") )) )    
        .def( bp::init< SireMol::EditMol const &, SireMol::ResID >(( bp::arg("editmol"), bp::arg("resid") )) )    
        .def( bp::init< SireMol::EditMol const &, QString const & >(( bp::arg("editmol"), bp::arg("resname") )) )    
        .def( 
            "add"
            , (void ( ::SireMol::EditRes::* )( ::QString const & ) )( &::SireMol::EditRes::add )
            , ( bp::arg("atm") ) )    
        .def( 
            "add"
            , (void ( ::SireMol::EditRes::* )( ::SireMol::Atom const & ) )( &::SireMol::EditRes::add )
            , ( bp::arg("atm") ) )    
        .def( 
            "add"
            , (void ( ::SireMol::EditRes::* )( ::QString const &,::SireMol::CutGroupNum ) )( &::SireMol::EditRes::add )
            , ( bp::arg("atm"), bp::arg("cgnum") ) )    
        .def( 
            "add"
            , (void ( ::SireMol::EditRes::* )( ::SireMol::Atom const &,::SireMol::CutGroupNum ) )( &::SireMol::EditRes::add )
            , ( bp::arg("atm"), bp::arg("cgnum") ) )    
        .def( 
            "add"
            , (void ( ::SireMol::EditRes::* )( ::SireMol::Bond const & ) )( &::SireMol::EditRes::add )
            , ( bp::arg("bond") ) )    
        .def( 
            "addAutoBonds"
            , (void ( ::SireMol::EditRes::* )(  ) )( &::SireMol::EditRes::addAutoBonds ) )    
        .def( 
            "addAutoBonds"
            , (void ( ::SireMol::EditRes::* )( ::SireMol::BondAddingFunction const & ) )( &::SireMol::EditRes::addAutoBonds )
            , ( bp::arg("bondfunc") ) )    
        .def( 
            "addBond"
            , &::SireMol::EditRes::addBond
            , ( bp::arg("atmnam0"), bp::arg("atmnam1") ) )    
        .def( 
            "angle"
            , &::SireMol::EditRes::angle
            , ( bp::arg("ang") ) )    
        .def( 
            "assertSameResidue"
            , (void ( ::SireMol::EditRes::* )( ::SireMol::AtomIndex const & ) const)( &::SireMol::EditRes::assertSameResidue )
            , ( bp::arg("atom") ) )    
        .def( 
            "assertSameResidue"
            , (void ( ::SireMol::EditRes::* )( ::QSet<SireMol::AtomIndex> const & ) const)( &::SireMol::EditRes::assertSameResidue )
            , ( bp::arg("atoms") ) )    
        .def( 
            "assertSameResidue"
            , (void ( ::SireMol::EditRes::* )( ::SireMol::CutGroupNum ) const)( &::SireMol::EditRes::assertSameResidue )
            , ( bp::arg("cgnum") ) )    
        .def( 
            "assertSameResidue"
            , (void ( ::SireMol::EditRes::* )( ::SireMol::CutGroupID ) const)( &::SireMol::EditRes::assertSameResidue )
            , ( bp::arg("cgid") ) )    
        .def( 
            "assertSameResidue"
            , (void ( ::SireMol::EditRes::* )( ::SireMol::Bond const & ) const)( &::SireMol::EditRes::assertSameResidue )
            , ( bp::arg("bond") ) )    
        .def( 
            "at"
            , (::SireMol::Atom ( ::SireMol::EditRes::* )( ::SireMol::AtomID ) const)( &::SireMol::EditRes::at )
            , ( bp::arg("i") ) )    
        .def( 
            "at"
            , (::SireMol::Atom ( ::SireMol::EditRes::* )( ::QString const & ) const)( &::SireMol::EditRes::at )
            , ( bp::arg("atomname") ) )    
        .def( 
            "at"
            , (::SireMol::Atom ( ::SireMol::EditRes::* )( ::SireMol::AtomIndex const & ) const)( &::SireMol::EditRes::at )
            , ( bp::arg("atom") ) )    
        .def( 
            "at"
            , (::SireMol::CutGroup ( ::SireMol::EditRes::* )( ::SireMol::CutGroupID ) const)( &::SireMol::EditRes::at )
            , ( bp::arg("cgid") ) )    
        .def( 
            "at"
            , (::SireMol::CutGroup ( ::SireMol::EditRes::* )( ::SireMol::CutGroupNum ) const)( &::SireMol::EditRes::at )
            , ( bp::arg("cgnum") ) )    
        .def( 
            "atom"
            , (::SireMol::Atom ( ::SireMol::EditRes::* )( ::SireMol::AtomID ) const)( &::SireMol::EditRes::atom )
            , ( bp::arg("i") ) )    
        .def( 
            "atom"
            , (::SireMol::Atom ( ::SireMol::EditRes::* )( ::QString const & ) const)( &::SireMol::EditRes::atom )
            , ( bp::arg("atomname") ) )    
        .def( 
            "atom"
            , (::SireMol::Atom ( ::SireMol::EditRes::* )( ::SireMol::AtomIndex const & ) const)( &::SireMol::EditRes::atom )
            , ( bp::arg("atom") ) )    
        .def( 
            "atomNames"
            , &::SireMol::EditRes::atomNames )    
        .def( 
            "atoms"
            , (::QVector<SireMol::Atom> ( ::SireMol::EditRes::* )(  ) const)( &::SireMol::EditRes::atoms ) )    
        .def( 
            "atoms"
            , (::QHash<SireMol::AtomID,SireMol::Atom> ( ::SireMol::EditRes::* )( ::QSet<SireMol::AtomID> const & ) const)( &::SireMol::EditRes::atoms )
            , ( bp::arg("atomids") ) )    
        .def( 
            "atoms"
            , (::QHash<QString,SireMol::Atom> ( ::SireMol::EditRes::* )( ::QSet<QString> const & ) const)( &::SireMol::EditRes::atoms )
            , ( bp::arg("atomnames") ) )    
        .def( 
            "atoms"
            , (::QHash<SireMol::AtomIndex,SireMol::Atom> ( ::SireMol::EditRes::* )( ::QSet<SireMol::AtomIndex> const & ) const)( &::SireMol::EditRes::atoms )
            , ( bp::arg("atms") ) )    
        .def( 
            "bond"
            , &::SireMol::EditRes::bond
            , ( bp::arg("bnd") ) )    
        .def( 
            "bondedResidues"
            , &::SireMol::EditRes::bondedResidues )    
        .def( 
            "connectivity"
            , &::SireMol::EditRes::connectivity )    
        .def( 
            "contains"
            , (bool ( ::SireMol::EditRes::* )( ::SireMol::CutGroupID ) const)( &::SireMol::EditRes::contains )
            , ( bp::arg("cgid") ) )    
        .def( 
            "contains"
            , (bool ( ::SireMol::EditRes::* )( ::SireMol::CutGroupNum ) const)( &::SireMol::EditRes::contains )
            , ( bp::arg("cgnum") ) )    
        .def( 
            "contains"
            , (bool ( ::SireMol::EditRes::* )( ::SireMol::AtomID ) const)( &::SireMol::EditRes::contains )
            , ( bp::arg("atomid") ) )    
        .def( 
            "contains"
            , (bool ( ::SireMol::EditRes::* )( ::QString const & ) const)( &::SireMol::EditRes::contains )
            , ( bp::arg("atomname") ) )    
        .def( 
            "contains"
            , (bool ( ::SireMol::EditRes::* )( ::SireMol::AtomIndex const & ) const)( &::SireMol::EditRes::contains )
            , ( bp::arg("atom") ) )    
        .def( 
            "contains"
            , (bool ( ::SireMol::EditRes::* )( ::SireMol::Bond const & ) const)( &::SireMol::EditRes::contains )
            , ( bp::arg("bond") ) )    
        .def( 
            "coordGroup"
            , (::SireVol::CoordGroup ( ::SireMol::EditRes::* )( ::SireMol::CutGroupID ) const)( &::SireMol::EditRes::coordGroup )
            , ( bp::arg("cgid") ) )    
        .def( 
            "coordGroup"
            , (::SireVol::CoordGroup ( ::SireMol::EditRes::* )( ::SireMol::CutGroupNum ) const)( &::SireMol::EditRes::coordGroup )
            , ( bp::arg("cgnum") ) )    
        .def( 
            "coordGroupsByID"
            , &::SireMol::EditRes::coordGroupsByID )    
        .def( 
            "coordGroupsByNum"
            , &::SireMol::EditRes::coordGroupsByNum )    
        .def( 
            "coordinates"
            , (::SireMaths::Vector ( ::SireMol::EditRes::* )( ::SireMol::AtomID ) const)( &::SireMol::EditRes::coordinates )
            , ( bp::arg("i") ) )    
        .def( 
            "coordinates"
            , (::SireMaths::Vector ( ::SireMol::EditRes::* )( ::QString const & ) const)( &::SireMol::EditRes::coordinates )
            , ( bp::arg("atomname") ) )    
        .def( 
            "coordinates"
            , (::SireMaths::Vector ( ::SireMol::EditRes::* )( ::SireMol::AtomIndex const & ) const)( &::SireMol::EditRes::coordinates )
            , ( bp::arg("atom") ) )    
        .def( 
            "coordinates"
            , (::QVector<SireMaths::Vector> ( ::SireMol::EditRes::* )(  ) const)( &::SireMol::EditRes::coordinates ) )    
        .def( 
            "coordinates"
            , (::QHash<SireMol::AtomID,SireMaths::Vector> ( ::SireMol::EditRes::* )( ::QSet<SireMol::AtomID> const & ) const)( &::SireMol::EditRes::coordinates )
            , ( bp::arg("atomids") ) )    
        .def( 
            "coordinates"
            , (::QHash<QString,SireMaths::Vector> ( ::SireMol::EditRes::* )( ::QSet<QString> const & ) const)( &::SireMol::EditRes::coordinates )
            , ( bp::arg("atomnames") ) )    
        .def( 
            "coordinates"
            , (::QHash<SireMol::AtomIndex,SireMaths::Vector> ( ::SireMol::EditRes::* )( ::QSet<SireMol::AtomIndex> const & ) const)( &::SireMol::EditRes::coordinates )
            , ( bp::arg("atms") ) )    
        .def( 
            "cutGroup"
            , (::SireMol::CutGroup ( ::SireMol::EditRes::* )( ::SireMol::CutGroupID ) const)( &::SireMol::EditRes::cutGroup )
            , ( bp::arg("cgid") ) )    
        .def( 
            "cutGroup"
            , (::SireMol::CutGroup ( ::SireMol::EditRes::* )( ::SireMol::CutGroupNum ) const)( &::SireMol::EditRes::cutGroup )
            , ( bp::arg("cgnum") ) )    
        .def( 
            "cutGroupsByID"
            , &::SireMol::EditRes::cutGroupsByID )    
        .def( 
            "cutGroupsByNum"
            , &::SireMol::EditRes::cutGroupsByNum )    
        .def( 
            "dihedral"
            , &::SireMol::EditRes::dihedral
            , ( bp::arg("dih") ) )    
        .def( 
            "getWeight"
            , (double ( ::SireMol::EditRes::* )( ::QStringList const &,::QStringList const &,::SireMol::WeightFunction const & ) const)( &::SireMol::EditRes::getWeight )
            , ( bp::arg("group0"), bp::arg("group1"), bp::arg("weightfunc") ) )    
        .def( 
            "getWeight"
            , (double ( ::SireMol::EditRes::* )( ::QSet<SireMol::AtomIndex> const &,::QSet<SireMol::AtomIndex> const &,::SireMol::WeightFunction const & ) const)( &::SireMol::EditRes::getWeight )
            , ( bp::arg("group0"), bp::arg("group1"), bp::arg("weightfunc") ) )    
        .def( 
            "improper"
            , &::SireMol::EditRes::improper
            , ( bp::arg("improper") ) )    
        .def( 
            "info"
            , &::SireMol::EditRes::info )    
        .def( 
            "isEmpty"
            , &::SireMol::EditRes::isEmpty )    
        .def( 
            "measure"
            , (double ( ::SireMol::EditRes::* )( ::SireMol::Bond const & ) const)( &::SireMol::EditRes::measure )
            , ( bp::arg("bnd") ) )    
        .def( 
            "measure"
            , (::SireUnits::Dimension::Angle ( ::SireMol::EditRes::* )( ::SireMol::Angle const & ) const)( &::SireMol::EditRes::measure )
            , ( bp::arg("ang") ) )    
        .def( 
            "measure"
            , (::SireUnits::Dimension::Angle ( ::SireMol::EditRes::* )( ::SireMol::Dihedral const & ) const)( &::SireMol::EditRes::measure )
            , ( bp::arg("dih") ) )    
        .def( 
            "measure"
            , (::SireUnits::Dimension::Angle ( ::SireMol::EditRes::* )( ::SireMol::Improper const & ) const)( &::SireMol::EditRes::measure )
            , ( bp::arg("improper") ) )    
        .def( 
            "molecule"
            , &::SireMol::EditRes::molecule )    
        .def( 
            "nAtoms"
            , (int ( ::SireMol::EditRes::* )(  ) const)( &::SireMol::EditRes::nAtoms ) )    
        .def( 
            "nAtoms"
            , (int ( ::SireMol::EditRes::* )( ::SireMol::CutGroupID ) const)( &::SireMol::EditRes::nAtoms )
            , ( bp::arg("cgid") ) )    
        .def( 
            "nAtoms"
            , (int ( ::SireMol::EditRes::* )( ::SireMol::CutGroupNum ) const)( &::SireMol::EditRes::nAtoms )
            , ( bp::arg("cgnum") ) )    
        .def( 
            "nBonds"
            , &::SireMol::EditRes::nBonds )    
        .def( 
            "nCutGroups"
            , &::SireMol::EditRes::nCutGroups )    
        .def( 
            "nInterBonds"
            , &::SireMol::EditRes::nInterBonds )    
        .def( 
            "nIntraBonds"
            , &::SireMol::EditRes::nIntraBonds )    
        .def( 
            "name"
            , &::SireMol::EditRes::name )    
        .def( 
            "number"
            , &::SireMol::EditRes::number )    
        .def( bp::self != bp::self )    
        .def( bp::self == bp::self )    
        .def( 
            "__getitem__"
            , (::SireMol::Atom ( ::SireMol::EditRes::* )( ::SireMol::AtomID ) const)( &::SireMol::EditRes::operator[] )
            , ( bp::arg("atomid") ) )    
        .def( 
            "__getitem__"
            , (::SireMol::Atom ( ::SireMol::EditRes::* )( ::QString const & ) const)( &::SireMol::EditRes::operator[] )
            , ( bp::arg("nam") ) )    
        .def( 
            "__getitem__"
            , (::SireMol::Atom ( ::SireMol::EditRes::* )( ::SireMol::AtomIndex const & ) const)( &::SireMol::EditRes::operator[] )
            , ( bp::arg("atom") ) )    
        .def( 
            "__getitem__"
            , (::SireMol::CutGroup ( ::SireMol::EditRes::* )( ::SireMol::CutGroupID ) const)( &::SireMol::EditRes::operator[] )
            , ( bp::arg("cgid") ) )    
        .def( 
            "__getitem__"
            , (::SireMol::CutGroup ( ::SireMol::EditRes::* )( ::SireMol::CutGroupNum ) const)( &::SireMol::EditRes::operator[] )
            , ( bp::arg("cgnum") ) )    
        .def( 
            "remove"
            , (void ( ::SireMol::EditRes::* )( ::QString const & ) )( &::SireMol::EditRes::remove )
            , ( bp::arg("atm") ) )    
        .def( 
            "remove"
            , (void ( ::SireMol::EditRes::* )( ::SireMol::AtomIndex const & ) )( &::SireMol::EditRes::remove )
            , ( bp::arg("atm") ) )    
        .def( 
            "remove"
            , (void ( ::SireMol::EditRes::* )( ::SireMol::Bond const & ) )( &::SireMol::EditRes::remove )
            , ( bp::arg("bond") ) )    
        .def( 
            "removeAllBonds"
            , (void ( ::SireMol::EditRes::* )( ::QString const & ) )( &::SireMol::EditRes::removeAllBonds )
            , ( bp::arg("atomname") ) )    
        .def( 
            "removeAllBonds"
            , (void ( ::SireMol::EditRes::* )(  ) )( &::SireMol::EditRes::removeAllBonds ) )    
        .def( 
            "removeBond"
            , &::SireMol::EditRes::removeBond
            , ( bp::arg("atmnam0"), bp::arg("atmnam1") ) )    
        .def( 
            "resName"
            , &::SireMol::EditRes::resName )    
        .def( 
            "resNum"
            , &::SireMol::EditRes::resNum )    
        .def( 
            "residuesBondedTo"
            , (::QHash<SireMol::ResNum,SireMol::EditRes> ( ::SireMol::EditRes::* )( ::SireMol::AtomID ) const)( &::SireMol::EditRes::residuesBondedTo )
            , ( bp::arg("atom") ) )    
        .def( 
            "residuesBondedTo"
            , (::QHash<SireMol::ResNum,SireMol::EditRes> ( ::SireMol::EditRes::* )( ::QString const & ) const)( &::SireMol::EditRes::residuesBondedTo )
            , ( bp::arg("atomname") ) )    
        .def( 
            "residuesBondedTo"
            , (::QHash<SireMol::ResNum,SireMol::EditRes> ( ::SireMol::EditRes::* )( ::SireMol::AtomIndex const & ) const)( &::SireMol::EditRes::residuesBondedTo )
            , ( bp::arg("atom") ) )    
        .def( 
            "setName"
            , &::SireMol::EditRes::setName
            , ( bp::arg("name") ) )    
        .def( 
            "setNumber"
            , &::SireMol::EditRes::setNumber
            , ( bp::arg("newnum") ) )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMol::EditRes >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMol::EditRes >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
