#ifndef SIREN_CLASS_H
#define SIREN_CLASS_H
/********************************************\
  *
  *  Siren - C++ metaobject library
  *
  *  Copyright (C) 2011  Christopher Woods
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

#include "Siren/siren.h"

SIREN_BEGIN_HEADER

namespace Siren
{

    /** This class is used to provide type-information about
        each of the Siren::Object derived classes
        
        @author Christopher Woods
    */
    class SIREN_EXPORT Class : public Object
    {
        SIREN_CLASS(Class,Object,1)

    public:
        Class();
        
        Class(const Object &object);

        Class(const char *type_name);
        Class(const String &type_name);
        
        Class(const Class &other);
        
        ~Class();

        String toString() const;
        
        String name() const;
        
        int version() const;
        
        bool hasSuper() const;
        
        Class superClass() const;
        
        bool implements(const String &type_name) const;
        bool implements(const char* type_name) const;

        template<class T>
        bool implements() const;
        
        StringList interfaces() const;
        
        Obj create() const;
        
        bool canCreate() const;
        bool isConcrete() const;
        
        bool canCast(const String &classname) const;
        bool canCast(const char *classname) const;

        void assertCanCast(const String &classname) const;
        void assertCanCast(const char *classname) const;
        
        template<class T>
        bool canCast() const;
        
        template<class T>
        void assertCanCast() const;

        static StringList registeredTypes();

    protected:
        friend class detail::ClassRegistry;
        Class(const detail::ClassData *metadata);

        void copy_object(const Class &other);
        bool compare_object(const Class &other) const;

    private:
        void checkNotNull() const;

        /** Pointer to the metatype data for this class */
        const detail::ClassData *d;

    }; // end of class Class

} // end of namespace Siren

SIREN_EXPOSE_CLASS( Siren::Class )

SIREN_END_HEADER

#endif
