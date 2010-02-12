// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License


#include "boost/python.hpp"

#include "Key.pypp.hpp"

#include "Lock.pypp.hpp"

#include "Password.pypp.hpp"

#include "PasswordLock.pypp.hpp"

namespace bp = boost::python;

#include "SireSec_properties.h"

BOOST_PYTHON_MODULE(_Sec){
    register_Key_class();

    register_Lock_class();

    register_Password_class();

    register_SireSec_properties();

    register_PasswordLock_class();
}
