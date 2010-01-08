// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "CoordGroupArray.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/quickcopy.hpp"

#include "SireMaths/axisset.h"

#include "SireMaths/matrix.h"

#include "SireMaths/quaternion.h"

#include "SireMaths/rotate.h"

#include "Siren/errors.h"

#include "Siren/stream.h"

#include "coordgroup.h"

#include <QDebug>

#include "coordgroup.h"

#include "Siren/logger.h"

SireVol::CoordGroupArray __copy__(const SireVol::CoordGroupArray &other){ return SireVol::CoordGroupArray(other); }

#include "Siren/str.hpp"

void register_CoordGroupArray_class(){

    { //::SireVol::CoordGroupArray
        typedef bp::class_< SireVol::CoordGroupArray, bp::bases< Siren::Mutable, Siren::Object > > CoordGroupArray_exposer_t;
        CoordGroupArray_exposer_t CoordGroupArray_exposer = CoordGroupArray_exposer_t( "CoordGroupArray" );
        bp::scope CoordGroupArray_scope( CoordGroupArray_exposer );
        CoordGroupArray_exposer.def( bp::init< >() );
        CoordGroupArray_exposer.def( bp::init< SireVol::CoordGroup const & >(( bp::arg("cgroup") )) );
        CoordGroupArray_exposer.def( bp::init< QVector< QVector< SireMaths::Vector > > const & >(( bp::arg("points") )) );
        CoordGroupArray_exposer.def( bp::init< QVector< SireVol::CoordGroup > const & >(( bp::arg("cgroups") )) );
        CoordGroupArray_exposer.def( bp::init< SireVol::CoordGroupArray const & >(( bp::arg("other") )) );
        { //::SireVol::CoordGroupArray::aaBox
        
            typedef ::SireVol::AABox ( ::SireVol::CoordGroupArray::*aaBox_function_type )(  ) const;
            aaBox_function_type aaBox_function_value( &::SireVol::CoordGroupArray::aaBox );
            
            CoordGroupArray_exposer.def( 
                "aaBox"
                , aaBox_function_value );
        
        }
        { //::SireVol::CoordGroupArray::assertValidCoordGroup
        
            typedef void ( ::SireVol::CoordGroupArray::*assertValidCoordGroup_function_type )( ::quint32 ) const;
            assertValidCoordGroup_function_type assertValidCoordGroup_function_value( &::SireVol::CoordGroupArray::assertValidCoordGroup );
            
            CoordGroupArray_exposer.def( 
                "assertValidCoordGroup"
                , assertValidCoordGroup_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireVol::CoordGroupArray::assertValidCoordinate
        
            typedef void ( ::SireVol::CoordGroupArray::*assertValidCoordinate_function_type )( ::quint32 ) const;
            assertValidCoordinate_function_type assertValidCoordinate_function_value( &::SireVol::CoordGroupArray::assertValidCoordinate );
            
            CoordGroupArray_exposer.def( 
                "assertValidCoordinate"
                , assertValidCoordinate_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireVol::CoordGroupArray::assertValidIndex
        
            typedef void ( ::SireVol::CoordGroupArray::*assertValidIndex_function_type )( ::quint32 ) const;
            assertValidIndex_function_type assertValidIndex_function_value( &::SireVol::CoordGroupArray::assertValidIndex );
            
            CoordGroupArray_exposer.def( 
                "assertValidIndex"
                , assertValidIndex_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireVol::CoordGroupArray::at
        
            typedef ::SireVol::CoordGroup const & ( ::SireVol::CoordGroupArray::*at_function_type )( ::quint32 ) const;
            at_function_type at_function_value( &::SireVol::CoordGroupArray::at );
            
            CoordGroupArray_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireVol::CoordGroupArray::changeFrame
        
            typedef void ( ::SireVol::CoordGroupArray::*changeFrame_function_type )( ::SireMaths::AxisSet const &,::SireMaths::AxisSet const & ) ;
            changeFrame_function_type changeFrame_function_value( &::SireVol::CoordGroupArray::changeFrame );
            
            CoordGroupArray_exposer.def( 
                "changeFrame"
                , changeFrame_function_value
                , ( bp::arg("from_frame"), bp::arg("to_frame") ) );
        
        }
        { //::SireVol::CoordGroupArray::changeFrame
        
            typedef void ( ::SireVol::CoordGroupArray::*changeFrame_function_type )( ::quint32,::SireMaths::AxisSet const &,::SireMaths::AxisSet const & ) ;
            changeFrame_function_type changeFrame_function_value( &::SireVol::CoordGroupArray::changeFrame );
            
            CoordGroupArray_exposer.def( 
                "changeFrame"
                , changeFrame_function_value
                , ( bp::arg("i"), bp::arg("from_frame"), bp::arg("to_frame") ) );
        
        }
        { //::SireVol::CoordGroupArray::count
        
            typedef int ( ::SireVol::CoordGroupArray::*count_function_type )(  ) const;
            count_function_type count_function_value( &::SireVol::CoordGroupArray::count );
            
            CoordGroupArray_exposer.def( 
                "count"
                , count_function_value );
        
        }
        { //::SireVol::CoordGroupArray::hashCode
        
            typedef ::uint ( ::SireVol::CoordGroupArray::*hashCode_function_type )(  ) const;
            hashCode_function_type hashCode_function_value( &::SireVol::CoordGroupArray::hashCode );
            
            CoordGroupArray_exposer.def( 
                "hashCode"
                , hashCode_function_value );
        
        }
        { //::SireVol::CoordGroupArray::mapInto
        
            typedef void ( ::SireVol::CoordGroupArray::*mapInto_function_type )( ::SireMaths::AxisSet const & ) ;
            mapInto_function_type mapInto_function_value( &::SireVol::CoordGroupArray::mapInto );
            
            CoordGroupArray_exposer.def( 
                "mapInto"
                , mapInto_function_value
                , ( bp::arg("axes") ) );
        
        }
        { //::SireVol::CoordGroupArray::mapInto
        
            typedef void ( ::SireVol::CoordGroupArray::*mapInto_function_type )( ::quint32,::SireMaths::AxisSet const & ) ;
            mapInto_function_type mapInto_function_value( &::SireVol::CoordGroupArray::mapInto );
            
            CoordGroupArray_exposer.def( 
                "mapInto"
                , mapInto_function_value
                , ( bp::arg("i"), bp::arg("axes") ) );
        
        }
        { //::SireVol::CoordGroupArray::merge
        
            typedef ::SireVol::CoordGroup ( ::SireVol::CoordGroupArray::*merge_function_type )(  ) const;
            merge_function_type merge_function_value( &::SireVol::CoordGroupArray::merge );
            
            CoordGroupArray_exposer.def( 
                "merge"
                , merge_function_value );
        
        }
        { //::SireVol::CoordGroupArray::nCoordGroups
        
            typedef int ( ::SireVol::CoordGroupArray::*nCoordGroups_function_type )(  ) const;
            nCoordGroups_function_type nCoordGroups_function_value( &::SireVol::CoordGroupArray::nCoordGroups );
            
            CoordGroupArray_exposer.def( 
                "nCoordGroups"
                , nCoordGroups_function_value );
        
        }
        { //::SireVol::CoordGroupArray::nCoords
        
            typedef int ( ::SireVol::CoordGroupArray::*nCoords_function_type )(  ) const;
            nCoords_function_type nCoords_function_value( &::SireVol::CoordGroupArray::nCoords );
            
            CoordGroupArray_exposer.def( 
                "nCoords"
                , nCoords_function_value );
        
        }
        CoordGroupArray_exposer.def( bp::self != bp::self );
        { //::SireVol::CoordGroupArray::operator=
        
            typedef ::SireVol::CoordGroupArray & ( ::SireVol::CoordGroupArray::*assign_function_type )( ::SireVol::CoordGroupArray const & ) ;
            assign_function_type assign_function_value( &::SireVol::CoordGroupArray::operator= );
            
            CoordGroupArray_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        CoordGroupArray_exposer.def( bp::self == bp::self );
        { //::SireVol::CoordGroupArray::operator[]
        
            typedef ::SireVol::CoordGroup const & ( ::SireVol::CoordGroupArray::*__getitem___function_type )( ::quint32 ) const;
            __getitem___function_type __getitem___function_value( &::SireVol::CoordGroupArray::operator[] );
            
            CoordGroupArray_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireVol::CoordGroupArray::restoreState
        
            typedef void ( ::SireVol::CoordGroupArray::*restoreState_function_type )( ::Siren::Object const & ) ;
            restoreState_function_type restoreState_function_value( &::SireVol::CoordGroupArray::restoreState );
            
            CoordGroupArray_exposer.def( 
                "restoreState"
                , restoreState_function_value
                , ( bp::arg("object") ) );
        
        }
        { //::SireVol::CoordGroupArray::rotate
        
            typedef void ( ::SireVol::CoordGroupArray::*rotate_function_type )( ::SireMaths::Quaternion const &,::SireMaths::Vector const & ) ;
            rotate_function_type rotate_function_value( &::SireVol::CoordGroupArray::rotate );
            
            CoordGroupArray_exposer.def( 
                "rotate"
                , rotate_function_value
                , ( bp::arg("quat"), bp::arg("point") ) );
        
        }
        { //::SireVol::CoordGroupArray::rotate
        
            typedef void ( ::SireVol::CoordGroupArray::*rotate_function_type )( ::SireMaths::Matrix const &,::SireMaths::Vector const & ) ;
            rotate_function_type rotate_function_value( &::SireVol::CoordGroupArray::rotate );
            
            CoordGroupArray_exposer.def( 
                "rotate"
                , rotate_function_value
                , ( bp::arg("rotmat"), bp::arg("point") ) );
        
        }
        { //::SireVol::CoordGroupArray::rotate
        
            typedef void ( ::SireVol::CoordGroupArray::*rotate_function_type )( ::quint32,::SireMaths::Quaternion const &,::SireMaths::Vector const & ) ;
            rotate_function_type rotate_function_value( &::SireVol::CoordGroupArray::rotate );
            
            CoordGroupArray_exposer.def( 
                "rotate"
                , rotate_function_value
                , ( bp::arg("i"), bp::arg("quat"), bp::arg("point") ) );
        
        }
        { //::SireVol::CoordGroupArray::rotate
        
            typedef void ( ::SireVol::CoordGroupArray::*rotate_function_type )( ::quint32,::SireMaths::Matrix const &,::SireMaths::Vector const & ) ;
            rotate_function_type rotate_function_value( &::SireVol::CoordGroupArray::rotate );
            
            CoordGroupArray_exposer.def( 
                "rotate"
                , rotate_function_value
                , ( bp::arg("i"), bp::arg("rotmat"), bp::arg("point") ) );
        
        }
        { //::SireVol::CoordGroupArray::saveState
        
            typedef ::Siren::ObjRef ( ::SireVol::CoordGroupArray::*saveState_function_type )(  ) const;
            saveState_function_type saveState_function_value( &::SireVol::CoordGroupArray::saveState );
            
            CoordGroupArray_exposer.def( 
                "saveState"
                , saveState_function_value );
        
        }
        { //::SireVol::CoordGroupArray::size
        
            typedef int ( ::SireVol::CoordGroupArray::*size_function_type )(  ) const;
            size_function_type size_function_value( &::SireVol::CoordGroupArray::size );
            
            CoordGroupArray_exposer.def( 
                "size"
                , size_function_value );
        
        }
        { //::SireVol::CoordGroupArray::stream
        
            typedef void ( ::SireVol::CoordGroupArray::*stream_function_type )( ::Siren::Stream & ) ;
            stream_function_type stream_function_value( &::SireVol::CoordGroupArray::stream );
            
            CoordGroupArray_exposer.def( 
                "stream"
                , stream_function_value
                , ( bp::arg("s") ) );
        
        }
        { //::SireVol::CoordGroupArray::toString
        
            typedef ::QString ( ::SireVol::CoordGroupArray::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireVol::CoordGroupArray::toString );
            
            CoordGroupArray_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireVol::CoordGroupArray::translate
        
            typedef void ( ::SireVol::CoordGroupArray::*translate_function_type )( ::SireMaths::Vector const & ) ;
            translate_function_type translate_function_value( &::SireVol::CoordGroupArray::translate );
            
            CoordGroupArray_exposer.def( 
                "translate"
                , translate_function_value
                , ( bp::arg("delta") ) );
        
        }
        { //::SireVol::CoordGroupArray::translate
        
            typedef void ( ::SireVol::CoordGroupArray::*translate_function_type )( ::quint32,::SireMaths::Vector const & ) ;
            translate_function_type translate_function_value( &::SireVol::CoordGroupArray::translate );
            
            CoordGroupArray_exposer.def( 
                "translate"
                , translate_function_value
                , ( bp::arg("i"), bp::arg("delta") ) );
        
        }
        { //::SireVol::CoordGroupArray::typeName
        
            typedef ::QString ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireVol::CoordGroupArray::typeName );
            
            CoordGroupArray_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireVol::CoordGroupArray::update
        
            typedef void ( ::SireVol::CoordGroupArray::*update_function_type )( ::quint32,::SireVol::CoordGroup const & ) ;
            update_function_type update_function_value( &::SireVol::CoordGroupArray::update );
            
            CoordGroupArray_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("i"), bp::arg("cgroup") ) );
        
        }
        { //::SireVol::CoordGroupArray::update
        
            typedef void ( ::SireVol::CoordGroupArray::*update_function_type )( ::quint32,::QVector< SireMaths::Vector > const & ) ;
            update_function_type update_function_value( &::SireVol::CoordGroupArray::update );
            
            CoordGroupArray_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("i"), bp::arg("coords") ) );
        
        }
        { //::SireVol::CoordGroupArray::update
        
            typedef void ( ::SireVol::CoordGroupArray::*update_function_type )( ::quint32,::SireMaths::Vector const *,int ) ;
            update_function_type update_function_value( &::SireVol::CoordGroupArray::update );
            
            CoordGroupArray_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("i"), bp::arg("coords"), bp::arg("ncoords") ) );
        
        }
        CoordGroupArray_exposer.staticmethod( "typeName" );
        CoordGroupArray_exposer.def( "__copy__", &__copy__);
        CoordGroupArray_exposer.def( "__deepcopy__", &__copy__);
        CoordGroupArray_exposer.def( "clone", &__copy__);
        CoordGroupArray_exposer.def( "__str__", &__str__< ::SireVol::CoordGroupArray > );
        CoordGroupArray_exposer.def( "__repr__", &__str__< ::SireVol::CoordGroupArray > );
    }

}
