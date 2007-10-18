// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "IOBase.pypp.hpp"
#include "boost/python.hpp"
#include "sireio_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/editmol.h"
#include "SireDB/parameterdb.h"

namespace bp = boost::python;

using namespace SireMol;

const char* pvt_get_name(const SireIO::IOBase&){ return "SireIO::IOBase";}

void register_IOBase_class(){

    bp::class_< SireIO::IOBase, boost::noncopyable >( "IOBase", bp::no_init )    
        .def( 
            "read"
            , (::QList<SireMol::Molecule> ( ::SireIO::IOBase::* )( ::QString,::SireMol::CuttingFunction const & ) const)( &::SireIO::IOBase::read )
            , ( bp::arg("filename"), bp::arg("cutfunc")=ResidueCutting() ) )    
        .def( 
            "read"
            , (::QList<SireMol::Molecule> ( ::SireIO::IOBase::* )( ::QIODevice &,::SireMol::CuttingFunction const & ) const)( &::SireIO::IOBase::read )
            , ( bp::arg("dev"), bp::arg("cutfunc")=ResidueCutting() ) )    
        .def( 
            "write"
            , (void ( ::SireIO::IOBase::* )( ::QList<SireMol::Molecule> const &,::QString ) const)( &::SireIO::IOBase::write )
            , ( bp::arg("molecules"), bp::arg("filename") ) )    
        .def( 
            "write"
            , (void ( ::SireIO::IOBase::* )( ::QList<SireMol::EditMol> const &,::QString ) const)( &::SireIO::IOBase::write )
            , ( bp::arg("molecules"), bp::arg("filename") ) )    
        .def( 
            "write"
            , (void ( ::SireIO::IOBase::* )( ::QList<SireMol::Molecule> const &,::QIODevice & ) const)( &::SireIO::IOBase::write )
            , ( bp::arg("molecules"), bp::arg("dev") ) )    
        .def( 
            "write"
            , (void ( ::SireIO::IOBase::* )( ::QList<SireMol::EditMol> const &,::QIODevice & ) const)( &::SireIO::IOBase::write )
            , ( bp::arg("molecules"), bp::arg("dev") ) )    
        .def( 
            "write"
            , (void ( ::SireIO::IOBase::* )( ::SireMol::Molecule const &,::QString ) const)( &::SireIO::IOBase::write )
            , ( bp::arg("molecule"), bp::arg("filename") ) )    
        .def( 
            "write"
            , (void ( ::SireIO::IOBase::* )( ::SireMol::Molecule const &,::QIODevice & ) const)( &::SireIO::IOBase::write )
            , ( bp::arg("molecule"), bp::arg("dev") ) )    
        .def( 
            "write"
            , (void ( ::SireIO::IOBase::* )( ::SireMol::EditMol const &,::QString ) const)( &::SireIO::IOBase::write )
            , ( bp::arg("molecule"), bp::arg("filename") ) )    
        .def( 
            "write"
            , (void ( ::SireIO::IOBase::* )( ::SireMol::EditMol const &,::QIODevice & ) const)( &::SireIO::IOBase::write )
            , ( bp::arg("molecule"), bp::arg("dev") ) )    
        .def( "__str__", &pvt_get_name);

}
