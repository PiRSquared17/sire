// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "_Base_free_functions.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "findexe.h"

#include <QDir>

#include <QProcess>

#include <QRegExp>

#include <QStringList>

#include "findexe.h"

#include "increment.h"

#include "increment.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "arrayproperty.hpp"

#include "numberproperty.h"

#include "propertylist.h"

#include "stringproperty.h"

#include "tostring.h"

#include "propertylist.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "arrayproperty.hpp"

#include "numberproperty.h"

#include "propertylist.h"

#include "stringproperty.h"

#include "tostring.h"

#include "propertylist.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "arrayproperty.hpp"

#include "numberproperty.h"

#include "propertylist.h"

#include "stringproperty.h"

#include "tostring.h"

#include "propertylist.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "arrayproperty.hpp"

#include "numberproperty.h"

#include "propertylist.h"

#include "stringproperty.h"

#include "tostring.h"

#include "propertylist.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "arrayproperty.hpp"

#include "numberproperty.h"

#include "propertylist.h"

#include "stringproperty.h"

#include "tostring.h"

#include "propertylist.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "arrayproperty.hpp"

#include "numberproperty.h"

#include "propertylist.h"

#include "stringproperty.h"

#include "tostring.h"

#include "propertylist.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "arrayproperty.hpp"

#include "numberproperty.h"

#include "propertylist.h"

#include "stringproperty.h"

#include "tostring.h"

#include "propertylist.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "arrayproperty.hpp"

#include "numberproperty.h"

#include "propertylist.h"

#include "stringproperty.h"

#include "tostring.h"

#include "propertylist.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "arrayproperty.hpp"

#include "numberproperty.h"

#include "propertylist.h"

#include "stringproperty.h"

#include "tostring.h"

#include "propertylist.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "arrayproperty.hpp"

#include "numberproperty.h"

#include "propertylist.h"

#include "stringproperty.h"

#include "tostring.h"

#include "propertylist.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "arrayproperty.hpp"

#include "numberproperty.h"

#include "propertylist.h"

#include "stringproperty.h"

#include "tostring.h"

#include "propertylist.h"

void register_free_functions(){

    { //::SireBase::findExe
    
        typedef ::QFileInfo ( *findExe_function_type )( ::QString const & );
        findExe_function_type findExe_function_value( &::SireBase::findExe );
        
        bp::def( 
            "findExe"
            , findExe_function_value
            , ( bp::arg("exe") ) );
    
    }

    { //::SireBase::increment
    
        typedef ::QString ( *increment_function_type )( ::QString );
        increment_function_type increment_function_value( &::SireBase::increment );
        
        bp::def( 
            "increment"
            , increment_function_value
            , ( bp::arg("name") ) );
    
    }

    { //::SireBase::wrap
    
        typedef ::SireBase::PropertyPtr ( *wrap_function_type )( ::QVector< QString > const & );
        wrap_function_type wrap_function_value( &::SireBase::wrap );
        
        bp::def( 
            "wrap"
            , wrap_function_value
            , ( bp::arg("values") ) );
    
    }

    { //::SireBase::wrap
    
        typedef ::SireBase::PropertyPtr ( *wrap_function_type )( ::QList< QString > const & );
        wrap_function_type wrap_function_value( &::SireBase::wrap );
        
        bp::def( 
            "wrap"
            , wrap_function_value
            , ( bp::arg("values") ) );
    
    }

    { //::SireBase::wrap
    
        typedef ::SireBase::PropertyPtr ( *wrap_function_type )( ::QVector< double > const & );
        wrap_function_type wrap_function_value( &::SireBase::wrap );
        
        bp::def( 
            "wrap"
            , wrap_function_value
            , ( bp::arg("values") ) );
    
    }

    { //::SireBase::wrap
    
        typedef ::SireBase::PropertyPtr ( *wrap_function_type )( ::QVector< int > const & );
        wrap_function_type wrap_function_value( &::SireBase::wrap );
        
        bp::def( 
            "wrap"
            , wrap_function_value
            , ( bp::arg("values") ) );
    
    }

    { //::SireBase::wrap
    
        typedef ::SireBase::PropertyPtr ( *wrap_function_type )( ::QList< double > const & );
        wrap_function_type wrap_function_value( &::SireBase::wrap );
        
        bp::def( 
            "wrap"
            , wrap_function_value
            , ( bp::arg("values") ) );
    
    }

    { //::SireBase::wrap
    
        typedef ::SireBase::PropertyPtr ( *wrap_function_type )( ::QList< int > const & );
        wrap_function_type wrap_function_value( &::SireBase::wrap );
        
        bp::def( 
            "wrap"
            , wrap_function_value
            , ( bp::arg("values") ) );
    
    }

    { //::SireBase::wrap
    
        typedef ::SireBase::PropertyPtr ( *wrap_function_type )( int );
        wrap_function_type wrap_function_value( &::SireBase::wrap );
        
        bp::def( 
            "wrap"
            , wrap_function_value
            , ( bp::arg("value") ) );
    
    }

    { //::SireBase::wrap
    
        typedef ::SireBase::PropertyPtr ( *wrap_function_type )( double );
        wrap_function_type wrap_function_value( &::SireBase::wrap );
        
        bp::def( 
            "wrap"
            , wrap_function_value
            , ( bp::arg("value") ) );
    
    }

    { //::SireBase::wrap
    
        typedef ::SireBase::PropertyPtr ( *wrap_function_type )( ::QString const & );
        wrap_function_type wrap_function_value( &::SireBase::wrap );
        
        bp::def( 
            "wrap"
            , wrap_function_value
            , ( bp::arg("value") ) );
    
    }

    { //::SireBase::wrap
    
        typedef ::SireBase::PropertyPtr ( *wrap_function_type )( ::QList< SireBase::PropPtr< SireBase::Property > > const & );
        wrap_function_type wrap_function_value( &::SireBase::wrap );
        
        bp::def( 
            "wrap"
            , wrap_function_value
            , ( bp::arg("value") ) );
    
    }

    { //::SireBase::wrap
    
        typedef ::SireBase::PropertyPtr ( *wrap_function_type )( ::SireBase::Property const & );
        wrap_function_type wrap_function_value( &::SireBase::wrap );
        
        bp::def( 
            "wrap"
            , wrap_function_value
            , ( bp::arg("value") ) );
    
    }

}
