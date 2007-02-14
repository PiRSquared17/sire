/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 2 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with this program; if not, write to the Free Software
  *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  *
  *  For full details of the license please see the COPYING file
  *  that should have come with this distribution.
  *
  *  You can contact the authors via the developer's mailing list
  *  at http://siremol.org
  *
\*********************************************/

#ifndef PYWRAP_SIREPY_CONVERTDICT_HPP
#define PYWRAP_SIREPY_CONVERTDICT_HPP

#include "sireglobal.h"

#include <boost/python.hpp>
#include <boost/tuple/tuple.hpp>

using namespace boost::python;

SIRE_BEGIN_HEADER

namespace SirePy
{

/** This struct provides the from-Python conversion from a dict
    to a dict or hash-like container of type 'C' (e.g. QHash, QMap) */
template<class C, class key_type, class value_type>
struct convert_python_to_dict
{
    /** Constructor - register the conversion functions
        for this type */
    convert_python_to_dict()
    {
        boost::python::converter::registry::push_back(
            &convertible,
            &construct,
            type_id< C >());
    }

    /** Test whether or not it is possible to convert the PyObject
        to a QVector where all of the elements are of type 'T' */
    static void* convertible(PyObject* obj_ptr)
    {
        //is this a dict type?
        if ( PyDict_Check(obj_ptr) )
        {
            //check the tuple elements... - convert to a boost::tuple object
            boost::python::dict d( handle<>(borrowed(obj_ptr)) );
            
            //check the items in the dict (items is a list of 2-tuples key-value
            list items = d.items();
            
            int nitems = extract<int>(items.attr("__len__")())();
            
            for (int i=0; i<nitems; ++i)
            {
                tuple item = extract<tuple>(items[i])();
                
                if ( not (extract<key_type>(item[0]).check() and 
                          extract<value_type>(item[1]).check()) )
                {
                    //either the key of value is wrong
                    return 0;
                }
            }
            
            //the tuple is ok!
            return obj_ptr;
        }
        else
            //could not recognise the type...
            return 0;
    }

    /** Construct a container of type T from the PyObject pointed to
        by 'obj_ptr' */
    static void construct(
        PyObject* obj_ptr,
        boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        //convert the PyObject to a boost::python::dict
        dict d( handle<>(borrowed(obj_ptr)) );
        
        //locate the storage space for the result
        void* storage =
            ( (converter::rvalue_from_python_storage<C>*)data )->storage.bytes;

        //create the container
        new (storage) C();

        C *container = static_cast<C*>(storage);

        //add all of the elements from the dict - do this by converting
        //to a list and then extracting each item
        list items = d.items();
            
        int nitems = extract<int>(items.attr("__len__")())();
        
        for (int i=0; i<nitems; ++i)
        {
            tuple item = extract<tuple>(items[i])();
                
            container->insert( extract<key_type>(item[0])(),
                               extract<value_type>(item[1])() );
        }

        data->convertible = storage;
    }
};

}

SIRE_END_HEADER

#endif
