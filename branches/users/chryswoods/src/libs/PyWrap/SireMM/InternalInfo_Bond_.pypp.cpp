// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "InternalInfo_Bond_.pypp.hpp"
#include "boost/python.hpp"
#include "siremm_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/atom.h"

namespace bp = boost::python;

SireMM::detail::InternalInfo<SireMol::Bond> __copy__(const SireMM::detail::InternalInfo<SireMol::Bond> &other){ return SireMM::detail::InternalInfo<SireMol::Bond>(other); }

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireMM::detail::InternalInfo<SireMol::Bond>&){ return "SireMM::detail::InternalInfo<SireMol::Bond>";}

void register_InternalInfo_Bond__class(){

    bp::class_< SireMM::detail::InternalInfo<SireMol::Bond> >( "InternalInfo_Bond_" )    
        .def( bp::init< >() )    
        .def( bp::init< SireMol::MoleculeInfo const & >(( bp::arg("molinfo") )) )    
        .def( 
            "at"
            , (::SireMol::GroupIndexID ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )( ::SireMol::Bond const & ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::at )
            , ( bp::arg("internal") ) )    
        .def( 
            "begin"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Bond> ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )(  ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::begin ) )    
        .def( 
            "contains"
            , (bool ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )( ::SireMol::Bond const & ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::contains )
            , ( bp::arg("internal") ) )    
        .def( 
            "contains"
            , (bool ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )( ::SireMol::ResNum ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::contains )
            , ( bp::arg("resnum") ) )    
        .def( 
            "contains"
            , (bool ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )( ::SireMol::AtomIndex const & ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::contains )
            , ( bp::arg("atom") ) )    
        .def( 
            "contains"
            , (bool ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )( ::SireMM::BondResID const & ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::contains )
            , ( bp::arg("resid") ) )    
        .def( 
            "contains"
            , (bool ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )( ::SireMol::GroupID ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::contains )
            , ( bp::arg("id") ) )    
        .def( 
            "contains"
            , (bool ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )( ::SireMol::GroupIndexID const & ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::contains )
            , ( bp::arg("id") ) )    
        .def( 
            "count"
            , (::size_t ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )(  ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::count ) )    
        .def( 
            "end"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Bond> ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )(  ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::end ) )    
        .def( 
            "find"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Bond> ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )( ::SireMol::Bond const & ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::find )
            , ( bp::arg("internal") ) )    
        .def( 
            "getInternal"
            , (::SireMol::Bond ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )( ::SireMol::GroupIndexID const & ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::getInternal )
            , ( bp::arg("id") ) )    
        .def( 
            "getInternals"
            , (::QSet<SireMol::Bond> const ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )( ::QSet<SireMol::GroupIndexID> const & ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::getInternals )
            , ( bp::arg("ids") ) )    
        .def( 
            "groupIDs"
            , (::QSet<SireMol::GroupID> const ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )(  ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::groupIDs ) )    
        .def( 
            "groupIDs"
            , (::QSet<SireMol::GroupID> const ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )( ::SireMol::ResNum ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::groupIDs )
            , ( bp::arg("resnum") ) )    
        .def( 
            "groupSizes"
            , (::QHash<SireMol::GroupID, unsigned> ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )(  ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::groupSizes ) )    
        .def( 
            "index"
            , (::SireMol::GroupIndexID ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )( ::SireMol::Bond const & ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::index )
            , ( bp::arg("internal") ) )    
        .def( 
            "info"
            , (::SireMol::MoleculeInfo const & ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )(  ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::info )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "interGroupIDs"
            , (::QSet<SireMol::GroupID> const ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )(  ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::interGroupIDs ) )    
        .def( 
            "interGroupIDs"
            , (::QSet<SireMol::GroupID> const ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )( ::SireMol::ResNum ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::interGroupIDs )
            , ( bp::arg("resnum") ) )    
        .def( 
            "interInternals"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Bond> ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )(  ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::interInternals ) )    
        .def( 
            "interInternals"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Bond> ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )( ::SireMol::ResNum ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::interInternals )
            , ( bp::arg("resnum") ) )    
        .def( 
            "internals"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Bond> ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )(  ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::internals ) )    
        .def( 
            "internals"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Bond> ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )( ::SireMol::ResNum ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::internals )
            , ( bp::arg("resnum") ) )    
        .def( 
            "internals"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Bond> ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )( ::SireMol::GroupID ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::internals )
            , ( bp::arg("groupid") ) )    
        .def( 
            "intraGroupIDs"
            , (::QSet<SireMol::GroupID> const ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )(  ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::intraGroupIDs ) )    
        .def( 
            "intraGroupIDs"
            , (::QSet<SireMol::GroupID> const ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )( ::SireMol::ResNum ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::intraGroupIDs )
            , ( bp::arg("resnum") ) )    
        .def( 
            "intraInternals"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Bond> ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )(  ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::intraInternals ) )    
        .def( 
            "intraInternals"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Bond> ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )( ::SireMol::ResNum ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::intraInternals )
            , ( bp::arg("resnum") ) )    
        .def( 
            "isEmpty"
            , (bool ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )(  ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::isEmpty ) )    
        .def( 
            "nGroups"
            , (int ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )(  ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::nGroups ) )    
        .def( 
            "nInterInternals"
            , (int ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )(  ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::nInterInternals ) )    
        .def( 
            "nInterInternals"
            , (int ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )( ::SireMol::ResNum ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::nInterInternals )
            , ( bp::arg("resnum") ) )    
        .def( 
            "nInternals"
            , (int ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )(  ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::nInternals ) )    
        .def( 
            "nInternals"
            , (int ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )( ::SireMol::GroupID ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::nInternals )
            , ( bp::arg("group") ) )    
        .def( 
            "nInternals"
            , (int ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )( ::SireMol::ResNum ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::nInternals )
            , ( bp::arg("resnum") ) )    
        .def( 
            "nIntraInternals"
            , (int ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )(  ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::nIntraInternals ) )    
        .def( 
            "nIntraInternals"
            , (int ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )( ::SireMol::ResNum ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::nIntraInternals )
            , ( bp::arg("resnum") ) )    
        .def( 
            "__getitem__"
            , (::SireMol::GroupIndexID ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )( ::SireMol::Bond const & ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::operator[] )
            , ( bp::arg("internal") ) )    
        .def( 
            "size"
            , (::size_t ( ::SireMM::detail::InternalInfo<SireMol::Bond>::* )(  ) const)( &::SireMM::detail::InternalInfo<SireMol::Bond>::size ) )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMM::detail::InternalInfo<SireMol::Bond> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMM::detail::InternalInfo<SireMol::Bond> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
