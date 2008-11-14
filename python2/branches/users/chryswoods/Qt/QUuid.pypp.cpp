// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "boost/python.hpp"
#include "sireqt_headers.h"
#include "QUuid.pypp.hpp"

#include <QUuid>

namespace bp = boost::python;

using namespace boost::python;

const char* pvt_get_name(const QUuid &uid){ return qPrintable(uid.toString());}

void register_QUuid_class()
{
    class_<QUuid>( "QUuid" )
       .def( "__str__", &pvt_get_name )
    ;

}