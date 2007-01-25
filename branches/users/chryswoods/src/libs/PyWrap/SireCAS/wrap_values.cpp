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

/**
  * This file contains the boost::python wrapping of Values 
  *
  * @author Christopher Woods
  *
  */
  
#include <Python.h>
#include <boost/python.hpp>
#include <QString>

#include "SireCAS/symbol.h"
#include "SireCAS/symbols.h"
#include "SireCAS/symbolvalue.h"
#include "SireCAS/identities.h"
#include "SireCAS/symbolexpression.h"
#include "SireCAS/values.h"

#include "SirePy/converttuple.hpp"

using namespace boost::python;

namespace SireCAS
{

void  
SIRECAS_EXPORT
export_Values()
{
    void (Values::*wrap_add0)(const SymbolValue&)
                    = &Values::add;
    void (Values::*wrap_add1)(const SymbolValue&, const SymbolValue&)
                    = &Values::add;
    void (Values::*wrap_add2)(const SymbolValue&, const SymbolValue&, const SymbolValue&)
                    = &Values::add;
    void (Values::*wrap_add3)(const SymbolValue&, const SymbolValue&, const SymbolValue&,
                              const SymbolValue&)
                    = &Values::add;
    void (Values::*wrap_add4)(const SymbolValue&, const SymbolValue&, const SymbolValue&,
                              const SymbolValue&, const SymbolValue&)
                    = &Values::add;
    void (Values::*wrap_add5)(const SymbolValue&, const SymbolValue&, const SymbolValue&,
                              const SymbolValue&, const SymbolValue&, const SymbolValue&)
                    = &Values::add;
    void (Values::*wrap_add6)(const SymbolValue&, const SymbolValue&, const SymbolValue&,
                              const SymbolValue&, const SymbolValue&, const SymbolValue&,
                              const SymbolValue&)
                    = &Values::add;
    void (Values::*wrap_add7)(const SymbolValue&, const SymbolValue&, const SymbolValue&,
                              const SymbolValue&, const SymbolValue&, const SymbolValue&,
                              const SymbolValue&, const SymbolValue&)
                    = &Values::add;
    void (Values::*wrap_add8)(const SymbolValue&, const SymbolValue&, const SymbolValue&,
                              const SymbolValue&, const SymbolValue&, const SymbolValue&,
                              const SymbolValue&, const SymbolValue&, const SymbolValue&)
                    = &Values::add;
    void (Values::*wrap_add9)(const SymbolValue&, const SymbolValue&, const SymbolValue&,
                              const SymbolValue&, const SymbolValue&, const SymbolValue&,
                              const SymbolValue&, const SymbolValue&, const SymbolValue&,
                              const SymbolValue&)
                    = &Values::add;
    
    class_<Values>("Values", init<>())
        .def( init<const QList<SymbolValue>&>() )
        .def( init<const SymbolValue&>() )
        .def( init<const SymbolValue&, const SymbolValue&>() )
        .def( init<const SymbolValue&, const SymbolValue&, const SymbolValue&>() )
        .def( init<const SymbolValue&, const SymbolValue&, const SymbolValue&,
                   const SymbolValue&>() )
        .def( init<const SymbolValue&, const SymbolValue&, const SymbolValue&,
                   const SymbolValue&, const SymbolValue&>() )
        .def( init<const SymbolValue&, const SymbolValue&, const SymbolValue&,
                   const SymbolValue&, const SymbolValue&, const SymbolValue&>() )
        .def( init<const SymbolValue&, const SymbolValue&, const SymbolValue&,
                   const SymbolValue&, const SymbolValue&, const SymbolValue&,
                   const SymbolValue&>() )
        .def( init<const SymbolValue&, const SymbolValue&, const SymbolValue&,
                   const SymbolValue&, const SymbolValue&, const SymbolValue&,
                   const SymbolValue&, const SymbolValue&>() )
        .def( init<const SymbolValue&, const SymbolValue&, const SymbolValue&,
                   const SymbolValue&, const SymbolValue&, const SymbolValue&,
                   const SymbolValue&, const SymbolValue&, const SymbolValue&>() )
        .def( init<const SymbolValue&, const SymbolValue&, const SymbolValue&,
                   const SymbolValue&, const SymbolValue&, const SymbolValue&,
                   const SymbolValue&, const SymbolValue&, const SymbolValue&,
                   const SymbolValue&>() )
        .def( "set", &Values::set )
        .def( "value", &Values::value )
        .def( "add", wrap_add0 )
        .def( "add", wrap_add1 )
        .def( "add", wrap_add2 )
        .def( "add", wrap_add3 )
        .def( "add", wrap_add4 )
        .def( "add", wrap_add5 )
        .def( "add", wrap_add6 )
        .def( "add", wrap_add7 )
        .def( "add", wrap_add8 )
        .def( "add", wrap_add9 )
    ;

    //allow implicit conversion from tuple(SymbolValue) to QList<SymbolValue>
    SirePy::convert_python_tuple_to_container< QList<SymbolValue> >();
    
    //allow implicit conversion from QList<SymbolValue> to Values
    implicitly_convertible<QList<SymbolValue>,Values>();

}

}
