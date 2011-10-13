#ifndef SIREN_MUTABLE_H
#define SIREN_MUTABLE_H
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
    /** This class provides the "Mutable" interface. The
        mutable interface must be attached to any Siren
        object that can be edited (changed).
        
        @author Christopher Woods
    */
    template<class T>
    class SIREN_EXPORT Mutable
    {
    public:
        typedef typename T Editor;
        typedef typename Editor::View View;
        
        Editor edit() const
        {
            const View *view = static_cast<const View*>(this);
            return Editor(*view);
        }
        
        Editor quickEdit()
        {
            View *view = static_cast<View*>(this);
            Editor editor(view);
            
            //destroy the original
            view->operator=(View());
            
            return editor;
        }
        
        View& operator=(const Editor &editor)
        {
            View *view = static_cast<View*>(this);
            return view->operator=( editor.commit() );
        }
        
        bool operator==(const Editor &editor) const
        {
            const View *view = static_cast<const View*>(this);
            return view->operator==(editor.commit());
        }
        
        bool operator!=(const Editor &editor) const
        {
            return not Mutable<T>::operator==(editor);
        }
        
    protected:
        Mutable()
        {}
        
        ~Mutable()
        {}

        static const char* interfaceName()
        {
            return "Mutable";
        }
        
        static bool IS_INTERFACE()
        {
            return true;
        }

    }; // end of class Mutable

} // end of namespace Siren

SIREN_END_HEADER

#endif // ifndef SIREN_MUTABLE_H
