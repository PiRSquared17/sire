// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "_Units_free_functions.pypp.hpp"

namespace bp = boost::python;

#include "angle.h"

#include "units.h"

#include "angle.h"

#include "convert.h"

#include "dimensions.h"

#include <QDebug>

#include <QMap>

#include <QPair>

#include <QString>

#include <QStringList>

#include <boost/shared_ptr.hpp>

#include "convert.h"

#include "convert.h"

#include "dimensions.h"

#include <QDebug>

#include <QMap>

#include <QPair>

#include <QString>

#include <QStringList>

#include <boost/shared_ptr.hpp>

#include "convert.h"

#include "convert.h"

#include "dimensions.h"

#include <QDebug>

#include <QMap>

#include <QPair>

#include <QString>

#include <QStringList>

#include <boost/shared_ptr.hpp>

#include "convert.h"

#include "convert.h"

#include "dimensions.h"

#include <QDebug>

#include <QMap>

#include <QPair>

#include <QString>

#include <QStringList>

#include <boost/shared_ptr.hpp>

#include "convert.h"

#include "convert.h"

#include "dimensions.h"

#include <QDebug>

#include <QMap>

#include <QPair>

#include <QString>

#include <QStringList>

#include <boost/shared_ptr.hpp>

#include "convert.h"

#include "convert.h"

#include "dimensions.h"

#include <QDebug>

#include <QMap>

#include <QPair>

#include <QString>

#include <QStringList>

#include <boost/shared_ptr.hpp>

#include "convert.h"

#include "convert.h"

#include "dimensions.h"

#include <QDebug>

#include <QMap>

#include <QPair>

#include <QString>

#include <QStringList>

#include <boost/shared_ptr.hpp>

#include "convert.h"

#include "convert.h"

#include "dimensions.h"

#include <QDebug>

#include <QMap>

#include <QPair>

#include <QString>

#include <QStringList>

#include <boost/shared_ptr.hpp>

#include "convert.h"

#include "convert.h"

#include "dimensions.h"

#include <QDebug>

#include <QMap>

#include <QPair>

#include <QString>

#include <QStringList>

#include <boost/shared_ptr.hpp>

#include "convert.h"

void register_free_functions(){

    { //::SireUnits::acute
    
        typedef ::SireUnits::Dimension::Angle ( *acute_function_type )( ::SireUnits::Dimension::Angle );
        acute_function_type acute_function_value( &::SireUnits::acute );
        
        bp::def( 
            "acute"
            , acute_function_value
            , ( bp::arg("angle") ) );
    
    }

    { //::SireUnits::convert
    
        typedef double ( *convert_function_type )( double,::SireUnits::Dimension::Unit const & );
        convert_function_type convert_function_value( &::SireUnits::convert );
        
        bp::def( 
            "convert"
            , convert_function_value
            , ( bp::arg("val"), bp::arg("to_units") ) );
    
    }

    { //::SireUnits::convert
    
        typedef double ( *convert_function_type )( double,::SireUnits::Dimension::detail::Temperature const & );
        convert_function_type convert_function_value( &::SireUnits::convert );
        
        bp::def( 
            "convert"
            , convert_function_value
            , ( bp::arg("val"), bp::arg("to_units") ) );
    
    }

    { //::SireUnits::convert
    
        typedef double ( *convert_function_type )( double,::SireUnits::Dimension::Unit const &,::SireUnits::Dimension::Unit const & );
        convert_function_type convert_function_value( &::SireUnits::convert );
        
        bp::def( 
            "convert"
            , convert_function_value
            , ( bp::arg("val"), bp::arg("from_units"), bp::arg("to_units") ) );
    
    }

    { //::SireUnits::convert
    
        typedef double ( *convert_function_type )( double,::SireUnits::Dimension::Unit const &,::SireUnits::Dimension::detail::Temperature const & );
        convert_function_type convert_function_value( &::SireUnits::convert );
        
        bp::def( 
            "convert"
            , convert_function_value
            , ( bp::arg("val"), bp::arg("from_units"), bp::arg("to_units") ) );
    
    }

    { //::SireUnits::convert
    
        typedef double ( *convert_function_type )( double,::SireUnits::Dimension::detail::Temperature const &,::SireUnits::Dimension::detail::Temperature const & );
        convert_function_type convert_function_value( &::SireUnits::convert );
        
        bp::def( 
            "convert"
            , convert_function_value
            , ( bp::arg("val"), bp::arg("from_units"), bp::arg("to_units") ) );
    
    }

    { //::SireUnits::convertFrom
    
        typedef double ( *convertFrom_function_type )( double,::SireUnits::Dimension::Unit const & );
        convertFrom_function_type convertFrom_function_value( &::SireUnits::convertFrom );
        
        bp::def( 
            "convertFrom"
            , convertFrom_function_value
            , ( bp::arg("val"), bp::arg("from_units") ) );
    
    }

    { //::SireUnits::convertFrom
    
        typedef double ( *convertFrom_function_type )( double,::SireUnits::Dimension::detail::Temperature const & );
        convertFrom_function_type convertFrom_function_value( &::SireUnits::convertFrom );
        
        bp::def( 
            "convertFrom"
            , convertFrom_function_value
            , ( bp::arg("val"), bp::arg("from_units") ) );
    
    }

    { //::SireUnits::convertTo
    
        typedef double ( *convertTo_function_type )( double,::SireUnits::Dimension::Unit const & );
        convertTo_function_type convertTo_function_value( &::SireUnits::convertTo );
        
        bp::def( 
            "convertTo"
            , convertTo_function_value
            , ( bp::arg("val"), bp::arg("to_units") ) );
    
    }

    { //::SireUnits::convertTo
    
        typedef double ( *convertTo_function_type )( double,::SireUnits::Dimension::detail::Temperature const & );
        convertTo_function_type convertTo_function_value( &::SireUnits::convertTo );
        
        bp::def( 
            "convertTo"
            , convertTo_function_value
            , ( bp::arg("val"), bp::arg("to_units") ) );
    
    }

}
