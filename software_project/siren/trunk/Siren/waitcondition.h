#ifndef SIREN_WAITCONDITION_H
#define SIREN_WAITCONDITION_H
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
#include "Siren/block.h"

SIREN_BEGIN_HEADER

namespace Siren
{
    class Mutex;
    class for_ages;

    namespace detail
    {
        class WaitConditionData;
    }

    /** This is a simple wait condition. A wait condition allows 
        a thread to pause (sleep), until it is woken up by another
        thread, or until the end of for_ages has been signalled */
    class SIREN_EXPORT WaitCondition : public Block
    {
    public:
        WaitCondition();
        WaitCondition(const WaitCondition &other);
        
        ~WaitCondition();
        
        WaitCondition& operator=(const WaitCondition &other);
        
        static const char* typeName(){ return "Siren::WaitCondition";}
        
        void wait();
        bool wait(unsigned long time);
        
        void wait(Mutex *mutex);
        bool wait(Mutex *mutex, unsigned long time);
        
        void wait(ReadWriteLock *lock);
        bool wait(ReadWriteLock *lock, unsigned long time);
        
        void wakeOne();
        void wakeAll();

    protected:
        friend class BlockRef;
        WaitCondition(const exp_shared_ptr<detail::BlockData>::Type &ptr);

        static bool isOfType(const exp_shared_ptr<detail::BlockData>::Type &ptr);

    private:
        detail::WaitConditionData *d;
        
    }; // end of class WaitCondition

}

SIREN_EXPOSE_CLASS( Siren::WaitCondition )

SIREN_END_HEADER

#endif // ifndef SIREN_WAITCONDITION_H