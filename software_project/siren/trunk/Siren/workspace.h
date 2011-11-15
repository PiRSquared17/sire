#ifndef SIREN_WORKSPACE_H
#define SIREN_WORKSPACE_H

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
#include "Siren/readwritelock.h"

SIREN_BEGIN_HEADER

namespace Siren
{
    /** A Receipt is given in response to a message passed via
        a WorkSpace. The Receipt can be used to confirm whether
        or not a message was received, and get when it was read 
        
        @author Christopher Woods
    */
    class Receipt
    {
    public:
        Receipt()
        {}
        
        ~Receipt()
        {}
    };

    /** Receipts are given when a message is broadcast to more
        then one receiver */
    class Receipts
    {
    public:
        Receipts()
        {}
        
        ~Receipts()
        {}
    };

    RENAME WorkSpace to WSpace and create WorkSpace that holds
    a shared pointer to a WSpace and provides the same API. That 
    way I have proper sharing of a WSpace object, so no unexpected
    deletions, and I can return WorkSpace objects, and work with
    WorkSpace objects

    /** This is the base class of all WorkSpace types. A WorkSpace
        provides a shared area that can be used by a group of worker
        threads to collaborate on the processing of a WorkPacket.
        A WorkSpace class type must be created for each WorkPacket
        type that can be processed in parallel, and is specific
        for that type.
        
        WorkSpaces are non-copyable and explicitly shared. Access
        to a WorkSpace must be protected via the read/write lock 
        that is supplied with that WorkSpace.
        
        @author Christopher Woods
    */
    class SIREN_EXPORT WorkSpace : public noncopyable
    {
    public:
        WorkSpace();
        WorkSpace(int worker_id);

        virtual ~WorkSpace();
        
        virtual const char* what() const=0;
        
        void flush();
        void setAutoFlush(bool on);
        bool isAutoFlush();
        
        virtual int arraySize(const String &key);
        
        virtual void save(const String &key, const Object &object);
        virtual Obj load(const String &key);
        
        virtual void save(const String &key, int idx, const Object &object);
        virtual Obj load(const String &key, int idx);
        
        virtual void save(const String &key, int idx, const List<Obj>::Type &objects);
                  
        virtual List<Obj>::Type load(const String &key, int idx, int count);
        
        virtual Obj receive(int worker_id, const String &key);
        virtual Receipt send(int worker_id, const String &key, const Object &object);
        virtual Receipts broadcast(const String &key, const Object &object);
        
        virtual bool hasObjectFrom(int worker_id, const String &key);
        virtual bool waitForObjectFrom(int worker_id, const String &key, int ms);
        
        template<class T>
        bool isA() const;
        
        template<class T>
        T& asA();
        
        template<class T>
        const T& asA() const;
        
    protected:
        void createKey(const String &key);
        void createArrayKey(const String &key, int count);
        
        void receivedFrom(int worker_id, const String &key, const Obj &object);
        
        ReadWriteLock& lock();
        
    private:
        /** The lock used to protect access to the data
            in this workspace */
        ReadWriteLock lk;
        
        /** The set of all arrays of saved objects */
        Hash<String,Vector<Obj>::Type>::Type obj_arrays;
        
        /** The received objects from other worker threads */
        Hash<int,Hash<String,Obj>::Type>::Type received_objs;
        
    }; // end of class WorkSpace

    /** Simple synonym for the explicitly shared pointer to a WorkSpace */
    typedef exp_shared_ptr<WorkSpace>::Type WorkSpacePtr;

} // end of namespace Siren

SIREN_EXPOSE_CLASS( Siren::WorkSpace )

SIREN_END_HEADER

#endif // ifndef SIREN_WORKSPACE_H
